use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdout, Command, Stdio};
use std::time::Instant;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum AnchorMode {
    Syncmer,
    Minimizer,
}

#[derive(Debug)]
struct Config {
    reads: Vec<PathBuf>,
    out_dir: PathBuf,
    anchor_mode: AnchorMode,
    k: usize,
    s: usize,
    syncmer_pos: Option<usize>,
    window: usize,
    min_anchor_spacing: usize,
    min_anchor_coverage: u32,
    min_edge_coverage: u32,
    min_branch_ratio: f64,
    max_edges_per_state: usize,
    dedup_kmer: usize,
    containment_ratio: f64,
    min_unitig_len: usize,
    min_tip_len: usize,
    min_link_support: u32,
    read_junction_links: bool,
    bidirectional_links: bool,
    junction_rescue_support: u32,
    min_read_len: usize,
    max_reads: Option<usize>,
    genome_size: Option<u64>,
    asm_coverage: Option<f64>,
    min_overlap: Option<usize>,
    iterations: Option<usize>,
    threads: usize,
    run_minimap2: bool,
}

#[derive(Debug, Clone)]
struct AnchorNode {
    seq: String,
    coverage: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct EdgeKey {
    from: usize,
    to: usize,
}

#[derive(Debug, Clone)]
struct EdgeStats {
    coverage: u32,
    total_span: u64,
    min_span: usize,
    max_span: usize,
    suffix: String,
}

#[derive(Debug, Clone, Copy)]
struct AnchorHit {
    state: usize,
    pos: usize,
}

#[derive(Debug)]
struct ReadWalk {
    edges: Vec<EdgeKey>,
}

#[derive(Debug)]
struct Assembly {
    nodes: Vec<AnchorNode>,
    key_to_node: HashMap<String, usize>,
    edges: HashMap<EdgeKey, EdgeStats>,
    edge_junctions: HashMap<(EdgeKey, EdgeKey), u32>,
    read_walks: Vec<ReadWalk>,
    reads_seen: usize,
    bases_seen: u64,
    anchors_seen: u64,
}

#[derive(Debug, Clone)]
struct Unitig {
    id: usize,
    path_states: Vec<usize>,
    path_edges: Vec<EdgeKey>,
    sequence: String,
    coverage: f64,
}

enum TextReader {
    Plain(BufReader<File>),
    Gzip {
        child: Child,
        reader: BufReader<ChildStdout>,
    },
}

impl Read for TextReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            TextReader::Plain(reader) => reader.read(buf),
            TextReader::Gzip { reader, .. } => reader.read(buf),
        }
    }
}

impl BufRead for TextReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            TextReader::Plain(reader) => reader.fill_buf(),
            TextReader::Gzip { reader, .. } => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            TextReader::Plain(reader) => reader.consume(amt),
            TextReader::Gzip { reader, .. } => reader.consume(amt),
        }
    }
}

impl Drop for TextReader {
    fn drop(&mut self) {
        if let TextReader::Gzip { child, .. } = self {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

fn main() {
    let started = Instant::now();
    let config = match Config::from_args(env::args().collect()) {
        Ok(config) => config,
        Err(message) => {
            eprintln!("{message}");
            print_usage();
            std::process::exit(2);
        }
    };

    if let Err(err) = run(config, started) {
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}

fn run(config: Config, started: Instant) -> io::Result<()> {
    validate_config(&config)?;
    fs::create_dir_all(&config.out_dir)?;

    let mut assembly = Assembly {
        nodes: Vec::new(),
        key_to_node: HashMap::new(),
        edges: HashMap::new(),
        edge_junctions: HashMap::new(),
        read_walks: Vec::new(),
        reads_seen: 0,
        bases_seen: 0,
        anchors_seen: 0,
    };

    let target_bases = match (config.genome_size, config.asm_coverage) {
        (Some(genome_size), Some(coverage)) if coverage > 0.0 => {
            Some((genome_size as f64 * coverage).round() as u64)
        }
        _ => None,
    };

    'inputs: for path in &config.reads {
        read_sequence_file(path, |name, seq| {
            if seq.len() < config.min_read_len {
                return Ok(true);
            }
            if let Some(max_reads) = config.max_reads {
                if assembly.reads_seen >= max_reads {
                    return Ok(false);
                }
            }
            if let Some(target) = target_bases {
                if assembly.bases_seen >= target {
                    return Ok(false);
                }
            }

            process_read(&config, &mut assembly, name, seq)?;
            Ok(true)
        })?;

        if let Some(max_reads) = config.max_reads {
            if assembly.reads_seen >= max_reads {
                break 'inputs;
            }
        }
        if let Some(target) = target_bases {
            if assembly.bases_seen >= target {
                break 'inputs;
            }
        }
    }

    let graph = build_filtered_graph(&config, &assembly);
    let compressed = compress_unitigs(&config, &assembly, &graph);
    let junctions = count_unitig_junctions(&assembly, &compressed.edge_to_unitig);
    let link_support = count_link_support(&assembly, &compressed.edge_to_placement);

    write_anchors(&config, &assembly, &graph)?;
    write_edges(&config, &assembly, &graph, &compressed.edge_to_unitig)?;
    write_unitigs(&config, &compressed.unitigs)?;
    write_gfa(
        &config,
        &compressed.unitigs,
        &compressed.state_starts,
        &compressed.state_ends,
        &link_support,
    )?;
    write_junctions(&config, &junctions)?;
    write_report(
        &config,
        &assembly,
        &graph,
        &compressed.unitigs,
        started.elapsed(),
    )?;

    if config.run_minimap2 {
        run_minimap2(&config)?;
    }

    Ok(())
}

fn validate_config(config: &Config) -> io::Result<()> {
    if config.reads.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "at least one --reads/-i/--pacbio-hifi input is required",
        ));
    }
    if config.k == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--k must be > 0",
        ));
    }
    if config.s == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--s must be > 0",
        ));
    }
    if config.anchor_mode == AnchorMode::Syncmer {
        if config.s > config.k {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--s must be <= --k in syncmer mode",
            ));
        }
        let t = config.syncmer_pos.unwrap_or((config.k - config.s) / 2);
        if t > config.k - config.s {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--syncmer-pos must be <= k - s",
            ));
        }
    }
    if config.window == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--window must be > 0",
        ));
    }
    if !(0.0..=1.0).contains(&config.min_branch_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-branch-ratio must be between 0 and 1",
        ));
    }
    if !(0.0..=1.0).contains(&config.containment_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--containment-ratio must be between 0 and 1",
        ));
    }
    Ok(())
}

impl Config {
    fn from_args(args: Vec<String>) -> Result<Self, String> {
        let mut config = Config {
            reads: Vec::new(),
            out_dir: PathBuf::from("result_graph"),
            anchor_mode: AnchorMode::Syncmer,
            k: 501,
            s: 31,
            syncmer_pos: None,
            window: 101,
            min_anchor_spacing: 0,
            min_anchor_coverage: 2,
            min_edge_coverage: 2,
            min_branch_ratio: 0.0,
            max_edges_per_state: 0,
            dedup_kmer: 0,
            containment_ratio: 0.0,
            min_unitig_len: 0,
            min_tip_len: 0,
            min_link_support: 0,
            read_junction_links: false,
            bidirectional_links: false,
            junction_rescue_support: 0,
            min_read_len: 0,
            max_reads: None,
            genome_size: None,
            asm_coverage: None,
            min_overlap: None,
            iterations: None,
            threads: 1,
            run_minimap2: false,
        };

        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "-i" | "--reads" | "--pacbio-hifi" => {
                    i += 1;
                    config
                        .reads
                        .push(PathBuf::from(take_arg(&args, i, "input path")?));
                }
                "-o" | "--out-dir" => {
                    i += 1;
                    config.out_dir = PathBuf::from(take_arg(&args, i, "output directory")?);
                }
                "--anchor" => {
                    i += 1;
                    config.anchor_mode = match take_arg(&args, i, "anchor mode")?.as_str() {
                        "syncmer" => AnchorMode::Syncmer,
                        "minimizer" => AnchorMode::Minimizer,
                        other => return Err(format!("unknown --anchor mode: {other}")),
                    };
                }
                "--k" | "-k" => {
                    i += 1;
                    config.k = parse_usize(take_arg(&args, i, "k")?, "--k")?;
                }
                "--s" => {
                    i += 1;
                    config.s = parse_usize(take_arg(&args, i, "s")?, "--s")?;
                }
                "--syncmer-pos" => {
                    i += 1;
                    config.syncmer_pos = Some(parse_usize(
                        take_arg(&args, i, "syncmer position")?,
                        "--syncmer-pos",
                    )?);
                }
                "--window" | "-w" => {
                    i += 1;
                    config.window = parse_usize(take_arg(&args, i, "window")?, "--window")?;
                }
                "--min-anchor-spacing" => {
                    i += 1;
                    config.min_anchor_spacing = parse_usize(
                        take_arg(&args, i, "minimum anchor spacing")?,
                        "--min-anchor-spacing",
                    )?;
                }
                "--min-anchor-coverage" => {
                    i += 1;
                    config.min_anchor_coverage = parse_u32(
                        take_arg(&args, i, "min anchor coverage")?,
                        "--min-anchor-coverage",
                    )?;
                }
                "--min-edge-coverage" => {
                    i += 1;
                    config.min_edge_coverage = parse_u32(
                        take_arg(&args, i, "min edge coverage")?,
                        "--min-edge-coverage",
                    )?;
                }
                "--min-branch-ratio" => {
                    i += 1;
                    config.min_branch_ratio = parse_f64(
                        take_arg(&args, i, "min branch ratio")?,
                        "--min-branch-ratio",
                    )?;
                }
                "--max-edges-per-state" => {
                    i += 1;
                    config.max_edges_per_state = parse_usize(
                        take_arg(&args, i, "max edges per state")?,
                        "--max-edges-per-state",
                    )?;
                }
                "--dedup-kmer" => {
                    i += 1;
                    config.dedup_kmer =
                        parse_usize(take_arg(&args, i, "dedup k-mer")?, "--dedup-kmer")?;
                }
                "--containment-ratio" => {
                    i += 1;
                    config.containment_ratio = parse_f64(
                        take_arg(&args, i, "containment ratio")?,
                        "--containment-ratio",
                    )?;
                }
                "--min-unitig-len" => {
                    i += 1;
                    config.min_unitig_len =
                        parse_usize(take_arg(&args, i, "min unitig length")?, "--min-unitig-len")?;
                }
                "--min-tip-len" => {
                    i += 1;
                    config.min_tip_len =
                        parse_usize(take_arg(&args, i, "min tip length")?, "--min-tip-len")?;
                }
                "--min-link-support" => {
                    i += 1;
                    config.min_link_support = parse_u32(
                        take_arg(&args, i, "min link support")?,
                        "--min-link-support",
                    )?;
                }
                "--read-junction-links" => {
                    config.read_junction_links = true;
                }
                "--bidirectional-links" => {
                    config.bidirectional_links = true;
                }
                "--junction-rescue-support" => {
                    i += 1;
                    config.junction_rescue_support = parse_u32(
                        take_arg(&args, i, "junction rescue support")?,
                        "--junction-rescue-support",
                    )?;
                }
                "--min-read-len" => {
                    i += 1;
                    config.min_read_len =
                        parse_usize(take_arg(&args, i, "min read length")?, "--min-read-len")?;
                }
                "--max-reads" => {
                    i += 1;
                    config.max_reads = Some(parse_usize(
                        take_arg(&args, i, "max reads")?,
                        "--max-reads",
                    )?);
                }
                "--genome-size" => {
                    i += 1;
                    config.genome_size = Some(parse_size(take_arg(&args, i, "genome size")?)?);
                }
                "--asm-coverage" => {
                    i += 1;
                    config.asm_coverage = Some(parse_f64(
                        take_arg(&args, i, "assembly coverage")?,
                        "--asm-coverage",
                    )?);
                }
                "--min-overlap" => {
                    i += 1;
                    config.min_overlap = Some(parse_usize(
                        take_arg(&args, i, "min overlap")?,
                        "--min-overlap",
                    )?);
                }
                "--iterations" => {
                    i += 1;
                    config.iterations = Some(parse_usize(
                        take_arg(&args, i, "iterations")?,
                        "--iterations",
                    )?);
                }
                "-t" | "--threads" => {
                    i += 1;
                    config.threads = parse_usize(take_arg(&args, i, "threads")?, "--threads")?;
                }
                "--minimap2" => {
                    config.run_minimap2 = true;
                }
                "-h" | "--help" => {
                    print_usage();
                    std::process::exit(0);
                }
                value if value.starts_with('-') => return Err(format!("unknown option: {value}")),
                value => {
                    config.reads.push(PathBuf::from(value));
                }
            }
            i += 1;
        }

        Ok(config)
    }
}

fn take_arg(args: &[String], index: usize, name: &str) -> Result<String, String> {
    args.get(index)
        .cloned()
        .ok_or_else(|| format!("missing value for {name}"))
}

fn parse_usize(value: String, flag: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

fn parse_u32(value: String, flag: &str) -> Result<u32, String> {
    value
        .parse::<u32>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

fn parse_f64(value: String, flag: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

fn parse_size(value: String) -> Result<u64, String> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Err("empty size".to_string());
    }
    let (number, multiplier) = match trimmed.as_bytes()[trimmed.len() - 1].to_ascii_lowercase() {
        b'k' => (&trimmed[..trimmed.len() - 1], 1_000_u64),
        b'm' => (&trimmed[..trimmed.len() - 1], 1_000_000_u64),
        b'g' => (&trimmed[..trimmed.len() - 1], 1_000_000_000_u64),
        _ => (trimmed, 1_u64),
    };
    let parsed = number
        .parse::<f64>()
        .map_err(|_| format!("invalid size: {value}"))?;
    Ok((parsed * multiplier as f64).round() as u64)
}

fn print_usage() {
    eprintln!(
        "Usage: simple_draft_asm --pacbio-hifi reads.fastq.gz -o result_graph [options]\n\
         \n\
         Inputs may be repeated with -i/--reads/--pacbio-hifi.\n\
         Key options:\n\
           --anchor syncmer|minimizer   default: syncmer\n\
           --k INT                      anchor k-mer length, default: 501\n\
           --s INT                      syncmer s-mer length, default: 31\n\
           --syncmer-pos INT            selected min s-mer offset, default: middle\n\
           --window INT                 minimizer window, default: 101\n\
           --min-anchor-spacing INT     thin selected anchors to this minimum read spacing\n\
           --genome-size 500k           enables --asm-coverage base limit\n\
           --asm-coverage FLOAT         target coverage limit\n\
           --min-edge-coverage INT      default: 2\n\
           --min-anchor-coverage INT    default: 2\n\
           --min-branch-ratio FLOAT     remove edges weak relative to local best edge\n\
           --max-edges-per-state INT    keep top N in/out edges per oriented anchor, 0 = unlimited\n\
           --dedup-kmer INT             k-mer size for post-compression containment pruning\n\
           --containment-ratio FLOAT    prune unitigs whose k-mer set is contained by a longer unitig\n\
           --min-unitig-len INT         drop compressed unitigs shorter than this length\n\
           --min-tip-len INT            drop terminal tips shorter than this length\n\
           --min-link-support INT       omit GFA links with lower read-walk support\n\
           --read-junction-links        emit additional read-walk-supported GFA links\n\
           --bidirectional-links        emit explicit reverse-complement GFA links\n\
           --junction-rescue-support INT keep early graph edges participating in supported raw junctions\n\
           --minimap2                   map reads to unitigs with external minimap2"
    );
}

fn read_sequence_file<F>(path: &Path, callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut reader = open_text(path)?;
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Ok(());
    }
    trim_line(&mut line);

    if line.starts_with('@') {
        read_fastq(reader, line, callback)
    } else if line.starts_with('>') {
        read_fasta(reader, line, callback)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} is not FASTA/FASTQ", path.display()),
        ))
    }
}

fn open_text(path: &Path) -> io::Result<TextReader> {
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let mut child = Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|err| {
                io::Error::new(
                    err.kind(),
                    format!("failed to spawn gzip for {}: {err}", path.display()),
                )
            })?;
        let stdout = child.stdout.take().ok_or_else(|| {
            io::Error::other(format!(
                "failed to capture gzip stdout for {}",
                path.display()
            ))
        })?;
        Ok(TextReader::Gzip {
            child,
            reader: BufReader::new(stdout),
        })
    } else {
        Ok(TextReader::Plain(BufReader::new(File::open(path)?)))
    }
}

fn read_fastq<F>(mut reader: TextReader, first_header: String, mut callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut header = first_header;
    loop {
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected FASTQ header, got {header}"),
            ));
        }
        let name = header[1..].split_whitespace().next().unwrap_or("");

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        if reader.read_line(&mut seq)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ sequence",
            ));
        }
        if reader.read_line(&mut plus)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ plus line",
            ));
        }
        if reader.read_line(&mut qual)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ quality",
            ));
        }
        trim_line(&mut seq);
        trim_line(&mut plus);
        trim_line(&mut qual);
        if !plus.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected FASTQ plus line for {name}"),
            ));
        }
        if !callback(name, &seq)? {
            break;
        }

        header.clear();
        if reader.read_line(&mut header)? == 0 {
            break;
        }
        trim_line(&mut header);
    }
    Ok(())
}

fn read_fasta<F>(mut reader: TextReader, first_header: String, mut callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut header = first_header;
    let mut seq = String::new();
    loop {
        let name = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        seq.clear();
        loop {
            let mut line = String::new();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                if !callback(&name, &seq)? {
                    return Ok(());
                }
                return Ok(());
            }
            trim_line(&mut line);
            if line.starts_with('>') {
                if !callback(&name, &seq)? {
                    return Ok(());
                }
                header = line;
                break;
            }
            seq.push_str(&line);
        }
    }
}

fn trim_line(line: &mut String) {
    while line.ends_with('\n') || line.ends_with('\r') {
        line.pop();
    }
}

fn process_read(
    config: &Config,
    assembly: &mut Assembly,
    _name: &str,
    seq: &str,
) -> io::Result<()> {
    assembly.reads_seen += 1;
    assembly.bases_seen += seq.len() as u64;

    let raw_hits = match config.anchor_mode {
        AnchorMode::Syncmer => {
            select_syncmer_hits(seq.as_bytes(), config.k, config.s, config.syncmer_pos)
        }
        AnchorMode::Minimizer => select_minimizer_hits(seq.as_bytes(), config.k, config.window),
    };
    let raw_hits = thin_anchor_hits(raw_hits, config.min_anchor_spacing);

    assembly.anchors_seen += raw_hits.len() as u64;

    let mut hits: Vec<AnchorHit> = Vec::with_capacity(raw_hits.len());
    let mut last_state: Option<usize> = None;
    for pos in raw_hits {
        let Some((canonical, is_forward)) = canonical_kmer(seq.as_bytes(), pos, config.k) else {
            continue;
        };
        let node_id = match assembly.key_to_node.get(&canonical) {
            Some(id) => *id,
            None => {
                let id = assembly.nodes.len();
                assembly.nodes.push(AnchorNode {
                    seq: canonical.clone(),
                    coverage: 0,
                });
                assembly.key_to_node.insert(canonical, id);
                id
            }
        };
        assembly.nodes[node_id].coverage += 1;
        let state = state_id(node_id, is_forward);
        if Some(state) != last_state {
            hits.push(AnchorHit { state, pos });
            last_state = Some(state);
        }
    }

    if hits.len() < 2 {
        assembly.read_walks.push(ReadWalk { edges: Vec::new() });
        return Ok(());
    }

    let seq_bytes = seq.as_bytes();
    let mut read_edges = Vec::with_capacity(hits.len().saturating_sub(1));
    for pair in hits.windows(2) {
        let left = pair[0];
        let right = pair[1];
        if right.pos <= left.pos {
            continue;
        }
        let suffix_start = left.pos.saturating_add(config.k).min(seq_bytes.len());
        let suffix_end = right.pos.saturating_add(config.k).min(seq_bytes.len());
        if suffix_end < suffix_start {
            continue;
        }
        let suffix =
            String::from_utf8_lossy(&seq_bytes[suffix_start..suffix_end]).to_ascii_uppercase();
        let span = right.pos.saturating_sub(left.pos);
        let edge_key = EdgeKey {
            from: left.state,
            to: right.state,
        };
        let entry = assembly.edges.entry(edge_key).or_insert_with(|| EdgeStats {
            coverage: 0,
            total_span: 0,
            min_span: span,
            max_span: span,
            suffix: suffix.clone(),
        });
        entry.coverage += 1;
        entry.total_span += span as u64;
        entry.min_span = entry.min_span.min(span);
        entry.max_span = entry.max_span.max(span);
        if suffix.len() > entry.suffix.len() {
            entry.suffix = suffix;
        }
        read_edges.push(edge_key);
    }

    for pair in read_edges.windows(2) {
        *assembly
            .edge_junctions
            .entry((pair[0], pair[1]))
            .or_insert(0) += 1;
    }

    assembly.read_walks.push(ReadWalk { edges: read_edges });
    Ok(())
}

fn select_syncmer_hits(seq: &[u8], k: usize, s: usize, syncmer_pos: Option<usize>) -> Vec<usize> {
    let target_offset = syncmer_pos.unwrap_or((k - s) / 2);
    let mut hits = Vec::new();

    for (segment_start, segment) in valid_segments(seq) {
        if segment.len() < k {
            continue;
        }
        let s_hashes = hash_smers(segment, s);
        let window = k - s + 1;
        let mut deque: VecDeque<usize> = VecDeque::new();

        for i in 0..s_hashes.len() {
            while let Some(&back) = deque.back() {
                if s_hashes[back] > s_hashes[i] || (s_hashes[back] == s_hashes[i] && back > i) {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(i);

            let k_start = (i + 1).saturating_sub(window);
            while let Some(&front) = deque.front() {
                if front < k_start {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if i + 1 >= window {
                let min_pos = *deque.front().expect("syncmer deque is not empty");
                if min_pos == k_start + target_offset {
                    hits.push(segment_start + k_start);
                }
            }
        }
    }

    hits
}

fn select_minimizer_hits(seq: &[u8], k: usize, window: usize) -> Vec<usize> {
    let mut hits = Vec::new();
    for (segment_start, segment) in valid_segments(seq) {
        if segment.len() < k {
            continue;
        }
        let hashes = hash_smers(segment, k);
        let mut deque: VecDeque<usize> = VecDeque::new();
        let mut emitted = None;

        for i in 0..hashes.len() {
            while let Some(&back) = deque.back() {
                if hashes[back] > hashes[i] || (hashes[back] == hashes[i] && back > i) {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(i);

            let window_start = (i + 1).saturating_sub(window);
            while let Some(&front) = deque.front() {
                if front < window_start {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if i + 1 >= window {
                let pos = *deque.front().expect("minimizer deque is not empty");
                if emitted != Some(pos) {
                    hits.push(segment_start + pos);
                    emitted = Some(pos);
                }
            }
        }
    }
    hits.sort_unstable();
    hits.dedup();
    hits
}

fn thin_anchor_hits(hits: Vec<usize>, min_spacing: usize) -> Vec<usize> {
    if min_spacing <= 1 || hits.len() < 2 {
        return hits;
    }

    let mut thinned = Vec::with_capacity(hits.len());
    let mut last_kept = None;
    for pos in hits {
        if last_kept.is_none_or(|last| pos >= last + min_spacing) {
            thinned.push(pos);
            last_kept = Some(pos);
        }
    }
    thinned
}

fn valid_segments(seq: &[u8]) -> Vec<(usize, &[u8])> {
    let mut segments = Vec::new();
    let mut start = None;
    for (i, &base) in seq.iter().enumerate() {
        if base_code(base).is_some() {
            if start.is_none() {
                start = Some(i);
            }
        } else if let Some(s) = start.take() {
            if i > s {
                segments.push((s, &seq[s..i]));
            }
        }
    }
    if let Some(s) = start {
        if seq.len() > s {
            segments.push((s, &seq[s..]));
        }
    }
    segments
}

fn hash_smers(segment: &[u8], s: usize) -> Vec<u64> {
    if s <= 32 {
        hash_smers_2bit(segment, s)
    } else {
        let mut hashes = Vec::with_capacity(segment.len().saturating_sub(s) + 1);
        for i in 0..=segment.len() - s {
            let forward = stable_hash_bytes(&segment[i..i + s]);
            let reverse = stable_hash_revcomp(&segment[i..i + s]);
            hashes.push(mix64(forward.min(reverse)));
        }
        hashes
    }
}

fn hash_smers_2bit(segment: &[u8], s: usize) -> Vec<u64> {
    let mask = if s == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * s)) - 1
    };
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mut hashes = Vec::with_capacity(segment.len().saturating_sub(s) + 1);

    for (i, &base) in segment.iter().enumerate() {
        let code = base_code(base).expect("valid segment contains only ACGT");
        forward = ((forward << 2) | code as u64) & mask;
        reverse = (reverse >> 2) | ((3 - code as u64) << (2 * (s - 1)));
        if i + 1 >= s {
            hashes.push(mix64(forward.min(reverse)));
        }
    }
    hashes
}

fn canonical_kmer(seq: &[u8], pos: usize, k: usize) -> Option<(String, bool)> {
    if pos + k > seq.len() {
        return None;
    }
    let slice = &seq[pos..pos + k];
    if !slice.iter().all(|&base| base_code(base).is_some()) {
        return None;
    }

    let mut forward = Vec::with_capacity(k);
    for &base in slice {
        forward.push(base.to_ascii_uppercase());
    }

    let mut reverse = Vec::with_capacity(k);
    for &base in slice.iter().rev() {
        reverse.push(complement(base.to_ascii_uppercase()));
    }

    if forward <= reverse {
        Some((String::from_utf8(forward).expect("ACGT is UTF-8"), true))
    } else {
        Some((String::from_utf8(reverse).expect("ACGT is UTF-8"), false))
    }
}

fn state_id(node_id: usize, is_forward: bool) -> usize {
    node_id * 2 + if is_forward { 0 } else { 1 }
}

fn node_id_from_state(state: usize) -> usize {
    state / 2
}

fn state_orientation(state: usize) -> char {
    if state.is_multiple_of(2) {
        '+'
    } else {
        '-'
    }
}

fn oriented_anchor_seq(node: &AnchorNode, state: usize) -> String {
    if state.is_multiple_of(2) {
        node.seq.clone()
    } else {
        revcomp_string(&node.seq)
    }
}

fn base_code(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn revcomp_string(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());
    for base in seq.bytes().rev() {
        out.push(complement(base) as char);
    }
    out
}

fn stable_hash_bytes(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for &byte in bytes {
        hash ^= byte.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(0x1000_0000_01b3);
    }
    hash
}

fn stable_hash_revcomp(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for &byte in bytes.iter().rev() {
        hash ^= complement(byte) as u64;
        hash = hash.wrapping_mul(0x1000_0000_01b3);
    }
    hash
}

fn mix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^ (x >> 33)
}

struct FilteredGraph {
    edges: Vec<(EdgeKey, EdgeStats)>,
    outgoing: HashMap<usize, Vec<EdgeKey>>,
    incoming: HashMap<usize, Vec<EdgeKey>>,
}

fn build_filtered_graph(config: &Config, assembly: &Assembly) -> FilteredGraph {
    let edge_junction_support = edge_junction_support(assembly);
    let mut candidates = Vec::new();
    let mut outgoing_rank: HashMap<usize, Vec<(EdgeKey, u32)>> = HashMap::new();
    let mut incoming_rank: HashMap<usize, Vec<(EdgeKey, u32)>> = HashMap::new();
    let mut max_outgoing: HashMap<usize, u32> = HashMap::new();
    let mut max_incoming: HashMap<usize, u32> = HashMap::new();

    for (key, stats) in &assembly.edges {
        let from_node = node_id_from_state(key.from);
        let to_node = node_id_from_state(key.to);
        let raw_junction_support = edge_junction_support.get(key).copied().unwrap_or(0);
        let rescue = config.junction_rescue_support > 0
            && raw_junction_support >= config.junction_rescue_support;
        if assembly.nodes[from_node].coverage < config.min_anchor_coverage {
            continue;
        }
        if assembly.nodes[to_node].coverage < config.min_anchor_coverage {
            continue;
        }
        if !rescue && stats.coverage < config.min_edge_coverage {
            continue;
        }

        candidates.push((*key, stats.clone()));
        let rank_support = stats.coverage.max(raw_junction_support);
        outgoing_rank
            .entry(key.from)
            .or_default()
            .push((*key, rank_support));
        incoming_rank
            .entry(key.to)
            .or_default()
            .push((*key, rank_support));
        max_outgoing
            .entry(key.from)
            .and_modify(|max_cov| *max_cov = (*max_cov).max(stats.coverage))
            .or_insert(stats.coverage);
        max_incoming
            .entry(key.to)
            .and_modify(|max_cov| *max_cov = (*max_cov).max(stats.coverage))
            .or_insert(stats.coverage);
    }

    let outgoing_keep = rank_keep_set(&mut outgoing_rank, config.max_edges_per_state);
    let incoming_keep = rank_keep_set(&mut incoming_rank, config.max_edges_per_state);
    let mut edges = Vec::new();
    let mut outgoing: HashMap<usize, Vec<EdgeKey>> = HashMap::new();
    let mut incoming: HashMap<usize, Vec<EdgeKey>> = HashMap::new();

    for (key, stats) in candidates {
        let raw_junction_support = edge_junction_support.get(&key).copied().unwrap_or(0);
        let rescue = config.junction_rescue_support > 0
            && raw_junction_support >= config.junction_rescue_support;
        if config.min_branch_ratio > 0.0 {
            let out_best = max_outgoing
                .get(&key.from)
                .copied()
                .unwrap_or(stats.coverage);
            let in_best = max_incoming.get(&key.to).copied().unwrap_or(stats.coverage);
            let out_ratio = stats.coverage as f64 / out_best.max(1) as f64;
            let in_ratio = stats.coverage as f64 / in_best.max(1) as f64;
            if !rescue
                && (out_ratio < config.min_branch_ratio || in_ratio < config.min_branch_ratio)
            {
                continue;
            }
        }
        if !rescue && !outgoing_keep.is_empty() && !outgoing_keep.contains(&key) {
            continue;
        }
        if !rescue && !incoming_keep.is_empty() && !incoming_keep.contains(&key) {
            continue;
        }
        edges.push((key, stats.clone()));
        outgoing.entry(key.from).or_default().push(key);
        incoming.entry(key.to).or_default().push(key);
    }

    edges.sort_by_key(|(key, _)| (key.from, key.to));
    for values in outgoing.values_mut() {
        values.sort_by_key(|key| key.to);
    }
    for values in incoming.values_mut() {
        values.sort_by_key(|key| key.from);
    }

    FilteredGraph {
        edges,
        outgoing,
        incoming,
    }
}

fn edge_junction_support(assembly: &Assembly) -> HashMap<EdgeKey, u32> {
    let mut support = HashMap::new();
    for (&(left, right), &count) in &assembly.edge_junctions {
        support
            .entry(left)
            .and_modify(|max_count: &mut u32| *max_count = (*max_count).max(count))
            .or_insert(count);
        support
            .entry(right)
            .and_modify(|max_count: &mut u32| *max_count = (*max_count).max(count))
            .or_insert(count);
    }
    support
}

fn rank_keep_set(
    ranks: &mut HashMap<usize, Vec<(EdgeKey, u32)>>,
    max_edges_per_state: usize,
) -> HashSet<EdgeKey> {
    let mut keep = HashSet::new();
    if max_edges_per_state == 0 {
        return keep;
    }

    for edges in ranks.values_mut() {
        edges.sort_by(|(left_key, left_cov), (right_key, right_cov)| {
            right_cov
                .cmp(left_cov)
                .then_with(|| left_key.from.cmp(&right_key.from))
                .then_with(|| left_key.to.cmp(&right_key.to))
        });
        for &(edge, _) in edges.iter().take(max_edges_per_state) {
            keep.insert(edge);
        }
    }

    keep
}

struct CompressedGraph {
    unitigs: Vec<Unitig>,
    edge_to_unitig: HashMap<EdgeKey, usize>,
    edge_to_placement: HashMap<EdgeKey, UnitigEnd>,
    state_starts: HashMap<usize, Vec<UnitigEnd>>,
    state_ends: HashMap<usize, Vec<UnitigEnd>>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct UnitigEnd {
    unitig_id: usize,
    orient: char,
}

fn compress_unitigs(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
) -> CompressedGraph {
    let mut visited: HashSet<EdgeKey> = HashSet::new();
    let mut directed_unitigs = Vec::new();

    let mut starts = Vec::new();
    for (edge, _) in &graph.edges {
        let indeg = graph.incoming.get(&edge.from).map_or(0, Vec::len);
        let outdeg = graph.outgoing.get(&edge.from).map_or(0, Vec::len);
        if indeg != 1 || outdeg != 1 {
            starts.push(*edge);
        }
    }

    for edge in starts {
        if visited.contains(&edge) {
            continue;
        }
        let unitig = trace_unitig(config, assembly, graph, edge, &mut visited);
        let unitig_id = directed_unitigs.len();
        directed_unitigs.push(Unitig {
            id: unitig_id,
            ..unitig
        });
    }

    for (edge, _) in &graph.edges {
        if visited.contains(edge) {
            continue;
        }
        let unitig = trace_unitig(config, assembly, graph, *edge, &mut visited);
        let unitig_id = directed_unitigs.len();
        directed_unitigs.push(Unitig {
            id: unitig_id,
            ..unitig
        });
    }

    let compressed = canonicalize_unitigs(directed_unitigs);
    let compressed = prune_redundant_unitigs(config, compressed);
    prune_short_tips(config, compressed)
}

fn canonicalize_unitigs(directed_unitigs: Vec<Unitig>) -> CompressedGraph {
    let mut path_index: HashMap<Vec<EdgeKey>, usize> = HashMap::new();
    for unitig in &directed_unitigs {
        path_index
            .entry(unitig.path_edges.clone())
            .or_insert(unitig.id);
    }

    let mut processed = vec![false; directed_unitigs.len()];
    let mut unitigs = Vec::new();
    let mut edge_to_unitig = HashMap::new();
    let mut edge_to_placement = HashMap::new();
    let mut state_starts: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();
    let mut state_ends: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();

    for old_id in 0..directed_unitigs.len() {
        if processed[old_id] {
            continue;
        }

        let old = &directed_unitigs[old_id];
        let rc_edges = reverse_complement_path(&old.path_edges);
        let rc_id = path_index.get(&rc_edges).copied();
        let use_old_forward = old.path_edges <= rc_edges;
        let segment_id = unitigs.len();
        let sequence = if use_old_forward {
            old.sequence.clone()
        } else {
            revcomp_string(&old.sequence)
        };
        let mut coverage_sum = old.coverage;
        let mut coverage_count = 1.0;

        processed[old_id] = true;
        let old_orient = if use_old_forward { '+' } else { '-' };
        register_unitig_projection(
            old,
            UnitigEnd {
                unitig_id: segment_id,
                orient: old_orient,
            },
            &mut edge_to_unitig,
            &mut edge_to_placement,
            &mut state_starts,
            &mut state_ends,
        );

        if let Some(pair_id) = rc_id {
            if pair_id != old_id {
                processed[pair_id] = true;
                let pair = &directed_unitigs[pair_id];
                coverage_sum += pair.coverage;
                coverage_count += 1.0;
                let pair_orient = if use_old_forward { '-' } else { '+' };
                register_unitig_projection(
                    pair,
                    UnitigEnd {
                        unitig_id: segment_id,
                        orient: pair_orient,
                    },
                    &mut edge_to_unitig,
                    &mut edge_to_placement,
                    &mut state_starts,
                    &mut state_ends,
                );
            }
        }

        unitigs.push(Unitig {
            id: segment_id,
            path_states: if use_old_forward {
                old.path_states.clone()
            } else {
                reverse_complement_states(&old.path_states)
            },
            path_edges: if use_old_forward {
                old.path_edges.clone()
            } else {
                rc_edges
            },
            sequence,
            coverage: coverage_sum / coverage_count,
        });
    }

    CompressedGraph {
        unitigs,
        edge_to_unitig,
        edge_to_placement,
        state_starts,
        state_ends,
    }
}

fn register_unitig_projection(
    unitig: &Unitig,
    placement: UnitigEnd,
    edge_to_unitig: &mut HashMap<EdgeKey, usize>,
    edge_to_placement: &mut HashMap<EdgeKey, UnitigEnd>,
    state_starts: &mut HashMap<usize, Vec<UnitigEnd>>,
    state_ends: &mut HashMap<usize, Vec<UnitigEnd>>,
) {
    for &edge in &unitig.path_edges {
        edge_to_unitig.insert(edge, placement.unitig_id);
        edge_to_placement.insert(edge, placement);
    }
    if let Some(&state) = unitig.path_states.first() {
        state_starts.entry(state).or_default().push(placement);
    }
    if let Some(&state) = unitig.path_states.last() {
        state_ends.entry(state).or_default().push(placement);
    }
}

fn reverse_complement_path(path_edges: &[EdgeKey]) -> Vec<EdgeKey> {
    path_edges
        .iter()
        .rev()
        .map(|edge| EdgeKey {
            from: reverse_state(edge.to),
            to: reverse_state(edge.from),
        })
        .collect()
}

fn reverse_complement_states(path_states: &[usize]) -> Vec<usize> {
    path_states
        .iter()
        .rev()
        .map(|&state| reverse_state(state))
        .collect()
}

fn reverse_state(state: usize) -> usize {
    state ^ 1
}

fn reverse_unitig_end(end: UnitigEnd) -> UnitigEnd {
    UnitigEnd {
        unitig_id: end.unitig_id,
        orient: if end.orient == '+' { '-' } else { '+' },
    }
}

fn prune_redundant_unitigs(config: &Config, graph: CompressedGraph) -> CompressedGraph {
    if config.min_unitig_len == 0
        && (config.dedup_kmer == 0 || config.containment_ratio <= 0.0 || graph.unitigs.len() < 2)
    {
        return graph;
    }

    let mut keep: Vec<bool> = graph
        .unitigs
        .iter()
        .map(|unitig| unitig.sequence.len() >= config.min_unitig_len)
        .collect();

    if config.dedup_kmer == 0 || config.containment_ratio <= 0.0 || graph.unitigs.len() < 2 {
        return rebuild_compressed_graph(graph, &keep);
    }

    let sketches: Vec<HashSet<u64>> = graph
        .unitigs
        .iter()
        .map(|unitig| sequence_kmer_set(&unitig.sequence, config.dedup_kmer))
        .collect();
    let mut order: Vec<usize> = (0..graph.unitigs.len()).collect();
    order.sort_by(|&left, &right| {
        graph.unitigs[right]
            .sequence
            .len()
            .cmp(&graph.unitigs[left].sequence.len())
            .then_with(|| left.cmp(&right))
    });

    for short_pos in 0..order.len() {
        let short_id = order[short_pos];
        if !keep[short_id] {
            continue;
        }
        if sketches[short_id].is_empty() {
            continue;
        }
        for &long_id in &order[..short_pos] {
            if !keep[long_id] || sketches[long_id].is_empty() {
                continue;
            }
            let shared = sketches[short_id]
                .iter()
                .filter(|hash| sketches[long_id].contains(hash))
                .count();
            let containment = shared as f64 / sketches[short_id].len() as f64;
            if containment >= config.containment_ratio {
                keep[short_id] = false;
                break;
            }
        }
    }

    rebuild_compressed_graph(graph, &keep)
}

fn prune_short_tips(config: &Config, graph: CompressedGraph) -> CompressedGraph {
    if config.min_tip_len == 0 || graph.unitigs.is_empty() {
        return graph;
    }

    let mut degrees = vec![0_usize; graph.unitigs.len()];
    let mut links = HashSet::new();
    for (state, end_ids) in &graph.state_ends {
        if let Some(start_ids) = graph.state_starts.get(state) {
            for &left in end_ids {
                for &right in start_ids {
                    if left.unitig_id == right.unitig_id && left.orient == right.orient {
                        continue;
                    }
                    if links.insert((left, right)) {
                        degrees[left.unitig_id] += 1;
                        degrees[right.unitig_id] += 1;
                    }
                }
            }
        }
    }

    let keep: Vec<bool> = graph
        .unitigs
        .iter()
        .map(|unitig| unitig.sequence.len() >= config.min_tip_len || degrees[unitig.id] > 1)
        .collect();
    rebuild_compressed_graph(graph, &keep)
}

fn sequence_kmer_set(sequence: &str, k: usize) -> HashSet<u64> {
    let seq = sequence.as_bytes();
    let mut hashes = HashSet::new();
    if k == 0 || seq.len() < k {
        return hashes;
    }
    for pos in 0..=seq.len() - k {
        if let Some(hash) = canonical_kmer_hash(seq, pos, k) {
            hashes.insert(hash);
        }
    }
    hashes
}

fn canonical_kmer_hash(seq: &[u8], pos: usize, k: usize) -> Option<u64> {
    if pos + k > seq.len() {
        return None;
    }
    let slice = &seq[pos..pos + k];
    if !slice.iter().all(|&base| base_code(base).is_some()) {
        return None;
    }
    let forward = stable_hash_bytes(slice);
    let reverse = stable_hash_revcomp(slice);
    Some(mix64(forward.min(reverse)))
}

fn rebuild_compressed_graph(graph: CompressedGraph, keep: &[bool]) -> CompressedGraph {
    let mut id_remap = HashMap::new();
    let mut unitigs = Vec::new();
    for unitig in graph.unitigs {
        if keep.get(unitig.id).copied().unwrap_or(false) {
            let new_id = unitigs.len();
            id_remap.insert(unitig.id, new_id);
            unitigs.push(Unitig {
                id: new_id,
                ..unitig
            });
        }
    }

    let edge_to_unitig = graph
        .edge_to_unitig
        .into_iter()
        .filter_map(|(edge, old_id)| id_remap.get(&old_id).map(|&new_id| (edge, new_id)))
        .collect();
    let edge_to_placement = graph
        .edge_to_placement
        .into_iter()
        .filter_map(|(edge, placement)| {
            id_remap.get(&placement.unitig_id).map(|&unitig_id| {
                (
                    edge,
                    UnitigEnd {
                        unitig_id,
                        orient: placement.orient,
                    },
                )
            })
        })
        .collect();
    let state_starts = remap_unitig_end_map(graph.state_starts, &id_remap);
    let state_ends = remap_unitig_end_map(graph.state_ends, &id_remap);

    CompressedGraph {
        unitigs,
        edge_to_unitig,
        edge_to_placement,
        state_starts,
        state_ends,
    }
}

fn remap_unitig_end_map(
    map: HashMap<usize, Vec<UnitigEnd>>,
    id_remap: &HashMap<usize, usize>,
) -> HashMap<usize, Vec<UnitigEnd>> {
    let mut remapped: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();
    let mut seen = HashSet::new();
    for (state, placements) in map {
        for placement in placements {
            let Some(&unitig_id) = id_remap.get(&placement.unitig_id) else {
                continue;
            };
            let new_placement = UnitigEnd {
                unitig_id,
                orient: placement.orient,
            };
            if seen.insert((state, new_placement)) {
                remapped.entry(state).or_default().push(new_placement);
            }
        }
    }
    remapped
}

fn trace_unitig(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
    first_edge: EdgeKey,
    visited: &mut HashSet<EdgeKey>,
) -> Unitig {
    let mut path_edges = Vec::new();
    let mut path_states = vec![first_edge.from];
    let mut sequence = oriented_anchor_seq(
        &assembly.nodes[node_id_from_state(first_edge.from)],
        first_edge.from,
    );
    let mut coverage_sum = 0_u64;
    let mut current = first_edge;

    loop {
        if visited.contains(&current) {
            break;
        }
        visited.insert(current);
        path_edges.push(current);
        path_states.push(current.to);

        if let Some(stats) = assembly.edges.get(&current) {
            sequence.push_str(&stats.suffix);
            coverage_sum += stats.coverage as u64;
        }

        let indeg = graph.incoming.get(&current.to).map_or(0, Vec::len);
        let out_edges = graph.outgoing.get(&current.to);
        let outdeg = out_edges.map_or(0, Vec::len);
        if indeg != 1 || outdeg != 1 {
            break;
        }
        let next = out_edges.expect("outdeg is 1")[0];
        if next == first_edge {
            break;
        }
        current = next;
    }

    let coverage = if path_edges.is_empty() {
        0.0
    } else {
        coverage_sum as f64 / path_edges.len() as f64
    };

    Unitig {
        id: usize::MAX,
        path_states,
        path_edges,
        sequence: if sequence.is_empty() {
            "N".repeat(config.k)
        } else {
            sequence
        },
        coverage,
    }
}

fn count_unitig_junctions(
    assembly: &Assembly,
    edge_to_unitig: &HashMap<EdgeKey, usize>,
) -> HashMap<(usize, usize), u32> {
    let mut counts = HashMap::new();
    for walk in &assembly.read_walks {
        let mut previous = None;
        for edge in &walk.edges {
            let Some(&unitig_id) = edge_to_unitig.get(edge) else {
                continue;
            };
            if previous != Some(unitig_id) {
                if let Some(prev_id) = previous {
                    *counts.entry((prev_id, unitig_id)).or_insert(0) += 1;
                }
                previous = Some(unitig_id);
            }
        }
    }
    counts
}

fn count_link_support(
    assembly: &Assembly,
    edge_to_placement: &HashMap<EdgeKey, UnitigEnd>,
) -> HashMap<(UnitigEnd, UnitigEnd), u32> {
    let mut counts = HashMap::new();
    for walk in &assembly.read_walks {
        let mut previous = None;
        for edge in &walk.edges {
            let Some(&placement) = edge_to_placement.get(edge) else {
                continue;
            };
            if previous != Some(placement) {
                if let Some(prev_placement) = previous {
                    *counts.entry((prev_placement, placement)).or_insert(0) += 1;
                }
                previous = Some(placement);
            }
        }
    }
    counts
}

fn write_anchors(config: &Config, assembly: &Assembly, graph: &FilteredGraph) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("anchors.tsv"))?;
    writeln!(
        out,
        "anchor_id\tcoverage\tforward_indegree\tforward_outdegree\treverse_indegree\treverse_outdegree\tsequence"
    )?;
    for (id, node) in assembly.nodes.iter().enumerate() {
        if node.coverage < config.min_anchor_coverage {
            continue;
        }
        let fwd = state_id(id, true);
        let rev = state_id(id, false);
        writeln!(
            out,
            "a{id}\t{}\t{}\t{}\t{}\t{}\t{}",
            node.coverage,
            graph.incoming.get(&fwd).map_or(0, Vec::len),
            graph.outgoing.get(&fwd).map_or(0, Vec::len),
            graph.incoming.get(&rev).map_or(0, Vec::len),
            graph.outgoing.get(&rev).map_or(0, Vec::len),
            node.seq
        )?;
    }
    Ok(())
}

fn write_edges(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
    edge_to_unitig: &HashMap<EdgeKey, usize>,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("edges.tsv"))?;
    let junction_support = edge_junction_support(assembly);
    writeln!(
        out,
        "edge_id\tfrom_anchor\tfrom_orient\tto_anchor\tto_orient\tcoverage\tmax_raw_junction_support\tmean_span\tmin_span\tmax_span\tunitig"
    )?;
    for (i, (edge, stats)) in graph.edges.iter().enumerate() {
        let mean_span = stats.total_span as f64 / stats.coverage.max(1) as f64;
        let raw_junction_support = junction_support.get(edge).copied().unwrap_or(0);
        let unitig = edge_to_unitig
            .get(edge)
            .map(|id| format!("utg{id}"))
            .unwrap_or_else(|| ".".to_string());
        writeln!(
            out,
            "e{i}\ta{}\t{}\ta{}\t{}\t{}\t{}\t{mean_span:.2}\t{}\t{}\t{}",
            node_id_from_state(edge.from),
            state_orientation(edge.from),
            node_id_from_state(edge.to),
            state_orientation(edge.to),
            stats.coverage,
            raw_junction_support,
            stats.min_span,
            stats.max_span,
            unitig
        )?;
    }
    Ok(())
}

fn write_unitigs(config: &Config, unitigs: &[Unitig]) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("unitigs.fasta"))?;
    for unitig in unitigs {
        writeln!(
            out,
            ">utg{} len={} edges={} cov={:.2}",
            unitig.id,
            unitig.sequence.len(),
            unitig.path_edges.len(),
            unitig.coverage
        )?;
        write_wrapped(&mut out, &unitig.sequence, 80)?;
    }
    Ok(())
}

fn write_wrapped(out: &mut File, sequence: &str, width: usize) -> io::Result<()> {
    for chunk in sequence.as_bytes().chunks(width) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn write_gfa(
    config: &Config,
    unitigs: &[Unitig],
    state_starts: &HashMap<usize, Vec<UnitigEnd>>,
    state_ends: &HashMap<usize, Vec<UnitigEnd>>,
    link_support: &HashMap<(UnitigEnd, UnitigEnd), u32>,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("graph.gfa"))?;
    writeln!(out, "H\tVN:Z:1.0")?;
    for unitig in unitigs {
        writeln!(
            out,
            "S\tutg{}\t{}\tLN:i:{}\tRC:f:{:.2}",
            unitig.id,
            unitig.sequence,
            unitig.sequence.len(),
            unitig.coverage
        )?;
    }

    let mut seen_links = HashSet::new();
    let mut links = Vec::new();
    for (state, end_ids) in state_ends {
        if let Some(start_ids) = state_starts.get(state) {
            for &left in end_ids {
                for &right in start_ids {
                    if left.unitig_id == right.unitig_id && left.orient == right.orient {
                        continue;
                    }
                    let support = link_support.get(&(left, right)).copied().unwrap_or(0);
                    if support < config.min_link_support {
                        continue;
                    }
                    add_gfa_link(
                        &mut links,
                        &mut seen_links,
                        left,
                        right,
                        support,
                        "terminal",
                        config.bidirectional_links,
                    );
                }
            }
        }
    }

    if config.read_junction_links {
        for (&(left, right), &support) in link_support {
            if left.unitig_id == right.unitig_id && left.orient == right.orient {
                continue;
            }
            if support < config.min_link_support {
                continue;
            }
            add_gfa_link(
                &mut links,
                &mut seen_links,
                left,
                right,
                support,
                "read",
                config.bidirectional_links,
            );
        }
    }

    links.sort_by_key(|(left, right, _support, source)| (*left, *right, *source));
    for (left, right, support, source) in links {
        if config.read_junction_links {
            writeln!(
                out,
                "L\tutg{}\t{}\tutg{}\t{}\t0M\tRC:i:{}\tJL:Z:{}",
                left.unitig_id, left.orient, right.unitig_id, right.orient, support, source
            )?;
        } else {
            writeln!(
                out,
                "L\tutg{}\t{}\tutg{}\t{}\t0M\tRC:i:{}",
                left.unitig_id, left.orient, right.unitig_id, right.orient, support
            )?;
        }
    }
    Ok(())
}

fn add_gfa_link(
    links: &mut Vec<(UnitigEnd, UnitigEnd, u32, &'static str)>,
    seen_links: &mut HashSet<(UnitigEnd, UnitigEnd)>,
    left: UnitigEnd,
    right: UnitigEnd,
    support: u32,
    source: &'static str,
    bidirectional: bool,
) {
    if seen_links.insert((left, right)) {
        links.push((left, right, support, source));
    }
    if bidirectional {
        let rc_left = reverse_unitig_end(right);
        let rc_right = reverse_unitig_end(left);
        if seen_links.insert((rc_left, rc_right)) {
            links.push((rc_left, rc_right, support, source));
        }
    }
}

fn write_junctions(config: &Config, junctions: &HashMap<(usize, usize), u32>) -> io::Result<()> {
    let mut out_totals: HashMap<usize, u32> = HashMap::new();
    for (&(from, _to), &count) in junctions {
        *out_totals.entry(from).or_insert(0) += count;
    }

    let mut rows: Vec<_> = junctions.iter().collect();
    rows.sort_by_key(|(&(from, to), _)| (from, to));

    let mut out = File::create(config.out_dir.join("junctions.tsv"))?;
    writeln!(out, "from_unitig\tto_unitig\tcount\tfrom_total\tfrequency")?;
    for (&(from, to), &count) in rows {
        let total = out_totals.get(&from).copied().unwrap_or(0);
        let freq = if total == 0 {
            0.0
        } else {
            count as f64 / total as f64
        };
        writeln!(out, "utg{from}\tutg{to}\t{count}\t{total}\t{freq:.6}")?;
    }
    Ok(())
}

fn write_report(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
    unitigs: &[Unitig],
    elapsed: std::time::Duration,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("report.txt"))?;
    writeln!(out, "simple_draft_asm report")?;
    writeln!(out, "elapsed_seconds\t{:.3}", elapsed.as_secs_f64())?;
    writeln!(out, "reads_seen\t{}", assembly.reads_seen)?;
    writeln!(out, "bases_seen\t{}", assembly.bases_seen)?;
    writeln!(out, "raw_anchor_hits\t{}", assembly.anchors_seen)?;
    writeln!(out, "unique_anchors\t{}", assembly.nodes.len())?;
    writeln!(out, "raw_edges\t{}", assembly.edges.len())?;
    writeln!(out, "raw_edge_junctions\t{}", assembly.edge_junctions.len())?;
    writeln!(out, "filtered_edges\t{}", graph.edges.len())?;
    writeln!(out, "unitigs\t{}", unitigs.len())?;
    writeln!(
        out,
        "unitig_bases\t{}",
        unitigs
            .iter()
            .map(|unitig| unitig.sequence.len())
            .sum::<usize>()
    )?;
    writeln!(out, "anchor_mode\t{:?}", config.anchor_mode)?;
    writeln!(out, "k\t{}", config.k)?;
    writeln!(out, "s\t{}", config.s)?;
    if config.anchor_mode == AnchorMode::Syncmer {
        writeln!(
            out,
            "syncmer_pos\t{}",
            config.syncmer_pos.unwrap_or((config.k - config.s) / 2)
        )?;
    } else {
        writeln!(out, "syncmer_pos\tNA")?;
    }
    writeln!(out, "window\t{}", config.window)?;
    writeln!(out, "min_anchor_spacing\t{}", config.min_anchor_spacing)?;
    writeln!(out, "min_anchor_coverage\t{}", config.min_anchor_coverage)?;
    writeln!(out, "min_edge_coverage\t{}", config.min_edge_coverage)?;
    writeln!(out, "min_branch_ratio\t{}", config.min_branch_ratio)?;
    writeln!(out, "max_edges_per_state\t{}", config.max_edges_per_state)?;
    writeln!(out, "dedup_kmer\t{}", config.dedup_kmer)?;
    writeln!(out, "containment_ratio\t{}", config.containment_ratio)?;
    writeln!(out, "min_unitig_len\t{}", config.min_unitig_len)?;
    writeln!(out, "min_tip_len\t{}", config.min_tip_len)?;
    writeln!(out, "min_link_support\t{}", config.min_link_support)?;
    writeln!(out, "read_junction_links\t{}", config.read_junction_links)?;
    writeln!(out, "bidirectional_links\t{}", config.bidirectional_links)?;
    writeln!(
        out,
        "junction_rescue_support\t{}",
        config.junction_rescue_support
    )?;
    if let Some(genome_size) = config.genome_size {
        writeln!(out, "genome_size\t{genome_size}")?;
    }
    if let Some(coverage) = config.asm_coverage {
        writeln!(out, "asm_coverage\t{coverage}")?;
    }
    if let Some(min_overlap) = config.min_overlap {
        writeln!(out, "min_overlap\t{min_overlap}")?;
    }
    if let Some(iterations) = config.iterations {
        writeln!(out, "iterations\t{iterations}")?;
    }
    writeln!(out, "threads\t{}", config.threads)?;
    Ok(())
}

fn run_minimap2(config: &Config) -> io::Result<()> {
    let unitigs = config.out_dir.join("unitigs.fasta");
    let paf = File::create(config.out_dir.join("read_to_unitigs.paf"))?;
    let mut command = Command::new("minimap2");
    command
        .arg("-x")
        .arg("map-hifi")
        .arg("-t")
        .arg(config.threads.to_string())
        .arg(unitigs);
    for read in &config.reads {
        command.arg(read);
    }
    let status = command.stdout(Stdio::from(paf)).status();
    match status {
        Ok(status) if status.success() => Ok(()),
        Ok(status) => Err(io::Error::other(format!(
            "minimap2 exited with status {status}"
        ))),
        Err(err) if err.kind() == io::ErrorKind::NotFound => {
            eprintln!("warning: --minimap2 requested, but minimap2 was not found in PATH");
            Ok(())
        }
        Err(err) => Err(err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_size_suffixes() {
        assert_eq!(parse_size("500k".to_string()).unwrap(), 500_000);
        assert_eq!(parse_size("1.5m".to_string()).unwrap(), 1_500_000);
    }

    #[test]
    fn reverse_complements() {
        assert_eq!(revcomp_string("ACGTTA"), "TAACGT");
    }

    #[test]
    fn canonicalizes_kmers() {
        let (canon, forward) = canonical_kmer(b"TTTACG", 0, 6).unwrap();
        assert_eq!(canon, "CGTAAA");
        assert!(!forward);
    }

    #[test]
    fn syncmer_selects_inside_valid_segment() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let hits = select_syncmer_hits(seq, 7, 3, Some(2));
        assert!(hits.iter().all(|&pos| pos + 7 <= seq.len()));
    }

    #[test]
    fn thins_anchor_hits_by_minimum_spacing() {
        assert_eq!(
            thin_anchor_hits(vec![2, 5, 13, 17, 28], 10),
            vec![2, 13, 28]
        );
    }
}
