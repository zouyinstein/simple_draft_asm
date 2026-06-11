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

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum OrganelleProfile {
    Plastid,
    Mito,
}

impl OrganelleProfile {
    fn as_str(self) -> &'static str {
        match self {
            OrganelleProfile::Plastid => "plastid",
            OrganelleProfile::Mito => "mito",
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum NumtInterference {
    Low,
    High,
}

impl NumtInterference {
    fn as_str(self) -> &'static str {
        match self {
            NumtInterference::Low => "low",
            NumtInterference::High => "high",
        }
    }
}

#[derive(Debug, Default)]
struct OverrideFlags {
    k: bool,
    s: bool,
    numt_interference: bool,
    min_anchor_coverage: bool,
    min_edge_coverage: bool,
    min_branch_ratio: bool,
    max_edges_per_state: bool,
    dedup_kmer: bool,
    containment_ratio: bool,
    min_tip_len: bool,
    min_link_support: bool,
    read_junction_links: bool,
    bidirectional_links: bool,
    rounds: bool,
}

#[derive(Debug, Clone)]
struct Config {
    reads: Vec<PathBuf>,
    out_dir: PathBuf,
    organelle: Option<OrganelleProfile>,
    numt_interference: NumtInterference,
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
    hifi_error_rate: f64,
    minimap_min_identity: f64,
    minimap_min_align_len: usize,
    paf_max_link_gap: isize,
    skeleton_gfa: Option<PathBuf>,
    skeleton_only: bool,
    skeleton_end_slop: usize,
    skeleton_min_link_support: u32,
    skeleton_min_link_ratio: f64,
    skeleton_rescue_gfa: Option<PathBuf>,
    skeleton_rescue_link_support: u32,
    rounds: usize,
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

    if config.skeleton_only {
        run_skeleton_linking(&config)?;
        return Ok(());
    }

    if config.rounds >= 2 {
        run_two_rounds(&config, started)?;
        return Ok(());
    }

    run_anchor_assembly(&config, started)
}

fn run_anchor_assembly(config: &Config, started: Instant) -> io::Result<()> {
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

    let full_config = full_graph_config(&config);
    let full_graph = build_filtered_graph(&full_config, &assembly);
    let full_compressed = compress_unitigs(&full_config, &assembly, &full_graph);
    let full_junctions = count_unitig_junctions(&assembly, &full_compressed.edge_to_unitig);
    let full_link_support = count_link_support(&assembly, &full_compressed.edge_to_placement);

    write_anchors(&config, &assembly, &graph)?;
    write_edges(&config, &assembly, &graph, &compressed.edge_to_unitig)?;
    write_unitigs(&config, &compressed.unitigs, "unitigs.fasta")?;
    write_unitigs(&config, &full_compressed.unitigs, "unitigs.full.fasta")?;
    write_gfa(
        &config,
        &compressed.unitigs,
        &compressed.state_starts,
        &compressed.state_ends,
        &link_support,
        "graph.gfa",
    )?;
    write_gfa(
        &config,
        &full_compressed.unitigs,
        &full_compressed.state_starts,
        &full_compressed.state_ends,
        &full_link_support,
        "graph.full.gfa",
    )?;
    write_junctions(&config, &junctions, "junctions.tsv")?;
    write_junctions(&config, &full_junctions, "junctions.full.tsv")?;
    write_report(
        &config,
        &assembly,
        &graph,
        &full_graph,
        &compressed.unitigs,
        &full_compressed.unitigs,
        started.elapsed(),
    )?;

    if config.run_minimap2 {
        run_minimap2(&config)?;
    }

    if config.skeleton_gfa.is_some() {
        run_skeleton_linking(&config)?;
    }

    Ok(())
}

fn run_two_rounds(config: &Config, started: Instant) -> io::Result<()> {
    fs::create_dir_all(&config.out_dir)?;

    let round1_dir = config.out_dir.join("round1_skeleton");
    let round1_readlinks_dir = config.out_dir.join("round1_readlinks");
    let round2_dir = config.out_dir.join("round2_skeleton");

    let mut round1 = config.clone();
    round1.out_dir = round1_dir.clone();
    round1.rounds = 1;
    round1.skeleton_gfa = None;
    round1.skeleton_rescue_gfa = None;
    round1.skeleton_only = false;
    round1.read_junction_links = false;
    round1.min_link_support = round1.min_link_support.min(10);
    round1.run_minimap2 = false;
    run_anchor_assembly(&round1, started)?;

    let mut round1_readlinks = round1.clone();
    round1_readlinks.out_dir = round1_readlinks_dir.clone();
    round1_readlinks.read_junction_links = true;
    run_anchor_assembly(&round1_readlinks, started)?;

    let mut round2 = config.clone();
    round2.out_dir = round2_dir.clone();
    round2.rounds = 1;
    round2.skeleton_only = true;
    round2.skeleton_gfa = Some(round1_dir.join("graph.gfa"));
    round2.skeleton_rescue_gfa = Some(round1_readlinks_dir.join("graph.gfa"));
    run_skeleton_linking(&round2)?;

    copy_final_two_round_outputs(config, &round1_dir, &round2_dir)?;
    write_two_round_report(
        config,
        &round1_dir,
        &round1_readlinks_dir,
        &round2_dir,
        started,
    )?;
    Ok(())
}

fn copy_final_two_round_outputs(
    config: &Config,
    round1_dir: &Path,
    round2_dir: &Path,
) -> io::Result<()> {
    copy_if_exists(
        round2_dir.join("skeleton.linked.gfa"),
        config.out_dir.join("graph.gfa"),
    )?;
    copy_if_exists(
        round2_dir.join("skeleton.links.tsv"),
        config.out_dir.join("links.tsv"),
    )?;
    copy_if_exists(
        round2_dir.join("skeleton.depth.tsv"),
        config.out_dir.join("depth.tsv"),
    )?;
    copy_if_exists(
        round2_dir.join("skeleton.report.txt"),
        config.out_dir.join("skeleton.report.txt"),
    )?;
    copy_if_exists(
        round2_dir.join("reads_to_skeleton.paf"),
        config.out_dir.join("reads_to_skeleton.paf"),
    )?;
    copy_if_exists(
        round1_dir.join("unitigs.fasta"),
        config.out_dir.join("unitigs.fasta"),
    )?;
    copy_if_exists(
        round1_dir.join("graph.gfa"),
        config.out_dir.join("graph.round1.gfa"),
    )?;
    copy_if_exists(
        round1_dir.join("graph.full.gfa"),
        config.out_dir.join("graph.round1.full.gfa"),
    )?;
    Ok(())
}

fn copy_if_exists(from: PathBuf, to: PathBuf) -> io::Result<()> {
    if from.exists() {
        fs::copy(&from, &to)?;
    }
    Ok(())
}

fn write_two_round_report(
    config: &Config,
    round1_dir: &Path,
    round1_readlinks_dir: &Path,
    round2_dir: &Path,
    started: Instant,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("report.txt"))?;
    writeln!(out, "simple_draft_asm two-round report")?;
    writeln!(
        out,
        "elapsed_seconds\t{:.3}",
        started.elapsed().as_secs_f64()
    )?;
    writeln!(
        out,
        "organelle\t{}",
        config
            .organelle
            .map(|o| o.as_str())
            .unwrap_or("unspecified")
    )?;
    writeln!(out, "rounds\t{}", config.rounds)?;
    writeln!(out, "round1_skeleton\t{}", round1_dir.display())?;
    writeln!(out, "round1_readlinks\t{}", round1_readlinks_dir.display())?;
    writeln!(out, "round2_skeleton\t{}", round2_dir.display())?;
    writeln!(out, "final_graph\tgraph.gfa")?;
    writeln!(out, "final_depth\tdepth.tsv")?;
    writeln!(out, "final_links\tlinks.tsv")?;
    writeln!(out, "round1_graph_copy\tgraph.round1.gfa")?;
    writeln!(out, "round1_full_graph_copy\tgraph.round1.full.gfa")?;
    writeln!(out, "k\t{}", config.k)?;
    writeln!(out, "s\t{}", config.s)?;
    writeln!(out, "min_anchor_coverage\t{}", config.min_anchor_coverage)?;
    writeln!(out, "min_edge_coverage\t{}", config.min_edge_coverage)?;
    writeln!(out, "min_branch_ratio\t{}", config.min_branch_ratio)?;
    writeln!(out, "max_edges_per_state\t{}", config.max_edges_per_state)?;
    writeln!(out, "dedup_kmer\t{}", config.dedup_kmer)?;
    writeln!(out, "containment_ratio\t{}", config.containment_ratio)?;
    writeln!(out, "min_tip_len\t{}", config.min_tip_len)?;
    writeln!(out, "min_link_support\t{}", config.min_link_support)?;
    writeln!(out, "skeleton_end_slop\t{}", config.skeleton_end_slop)?;
    writeln!(
        out,
        "skeleton_min_link_support\t{}",
        config.skeleton_min_link_support
    )?;
    writeln!(
        out,
        "skeleton_min_link_ratio\t{}",
        config.skeleton_min_link_ratio
    )?;
    writeln!(
        out,
        "skeleton_rescue_link_support\t{}",
        config.skeleton_rescue_link_support
    )?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    Ok(())
}

fn full_graph_config(config: &Config) -> Config {
    let mut full = config.clone();
    full.min_branch_ratio = 0.0;
    full.max_edges_per_state = 0;
    full.dedup_kmer = 0;
    full.containment_ratio = 0.0;
    full.min_unitig_len = 0;
    full.min_tip_len = 0;
    full
}

fn validate_config(config: &Config) -> io::Result<()> {
    if config.reads.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "at least one --reads/-i/--pacbio-hifi input is required",
        ));
    }
    if config.rounds == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--rounds must be >= 1",
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
    if !(0.0..1.0).contains(&config.hifi_error_rate) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--hifi-error-rate must be >= 0 and < 1",
        ));
    }
    if !(0.0..=1.0).contains(&config.minimap_min_identity) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--minimap-min-identity must be between 0 and 1",
        ));
    }
    if config.skeleton_only && config.skeleton_gfa.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--skeleton-only requires --skeleton-gfa",
        ));
    }
    if !(0.0..=1.0).contains(&config.skeleton_min_link_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--skeleton-min-link-ratio must be between 0 and 1",
        ));
    }
    Ok(())
}

impl Config {
    fn from_args(args: Vec<String>) -> Result<Self, String> {
        let mut config = Config {
            reads: Vec::new(),
            out_dir: PathBuf::from("result_graph"),
            organelle: None,
            numt_interference: NumtInterference::Low,
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
            hifi_error_rate: 0.003,
            minimap_min_identity: 0.95,
            minimap_min_align_len: 500,
            paf_max_link_gap: 2_000,
            skeleton_gfa: None,
            skeleton_only: false,
            skeleton_end_slop: 1_500,
            skeleton_min_link_support: 10,
            skeleton_min_link_ratio: 0.20,
            skeleton_rescue_gfa: None,
            skeleton_rescue_link_support: 20,
            rounds: 0,
            threads: 1,
            run_minimap2: false,
        };

        let mut overrides = OverrideFlags::default();
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
                "--organelle" => {
                    i += 1;
                    config.organelle = Some(match take_arg(&args, i, "organelle")?.as_str() {
                        "plastid" | "chloroplast" | "cp" => OrganelleProfile::Plastid,
                        "mito" | "mitochondria" | "mitochondrion" | "mt" => OrganelleProfile::Mito,
                        other => return Err(format!("unknown --organelle value: {other}")),
                    });
                }
                "--numt-interference" => {
                    i += 1;
                    config.numt_interference = match take_arg(&args, i, "NUMT interference")?
                        .as_str()
                    {
                        "low" | "small" | "minor" => NumtInterference::Low,
                        "high" | "large" | "strong" => NumtInterference::High,
                        other => return Err(format!("unknown --numt-interference value: {other}")),
                    };
                    overrides.numt_interference = true;
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
                    overrides.k = true;
                }
                "--s" => {
                    i += 1;
                    config.s = parse_usize(take_arg(&args, i, "s")?, "--s")?;
                    overrides.s = true;
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
                    overrides.min_anchor_coverage = true;
                }
                "--min-edge-coverage" => {
                    i += 1;
                    config.min_edge_coverage = parse_u32(
                        take_arg(&args, i, "min edge coverage")?,
                        "--min-edge-coverage",
                    )?;
                    overrides.min_edge_coverage = true;
                }
                "--min-branch-ratio" => {
                    i += 1;
                    config.min_branch_ratio = parse_f64(
                        take_arg(&args, i, "min branch ratio")?,
                        "--min-branch-ratio",
                    )?;
                    overrides.min_branch_ratio = true;
                }
                "--max-edges-per-state" => {
                    i += 1;
                    config.max_edges_per_state = parse_usize(
                        take_arg(&args, i, "max edges per state")?,
                        "--max-edges-per-state",
                    )?;
                    overrides.max_edges_per_state = true;
                }
                "--dedup-kmer" => {
                    i += 1;
                    config.dedup_kmer =
                        parse_usize(take_arg(&args, i, "dedup k-mer")?, "--dedup-kmer")?;
                    overrides.dedup_kmer = true;
                }
                "--containment-ratio" => {
                    i += 1;
                    config.containment_ratio = parse_f64(
                        take_arg(&args, i, "containment ratio")?,
                        "--containment-ratio",
                    )?;
                    overrides.containment_ratio = true;
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
                    overrides.min_tip_len = true;
                }
                "--min-link-support" => {
                    i += 1;
                    config.min_link_support = parse_u32(
                        take_arg(&args, i, "min link support")?,
                        "--min-link-support",
                    )?;
                    overrides.min_link_support = true;
                }
                "--read-junction-links" => {
                    config.read_junction_links = true;
                    overrides.read_junction_links = true;
                }
                "--bidirectional-links" => {
                    config.bidirectional_links = true;
                    overrides.bidirectional_links = true;
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
                "--hifi-error-rate" => {
                    i += 1;
                    config.hifi_error_rate =
                        parse_f64(take_arg(&args, i, "HiFi error rate")?, "--hifi-error-rate")?;
                }
                "--hifi-accuracy" => {
                    i += 1;
                    let accuracy =
                        parse_f64(take_arg(&args, i, "HiFi accuracy")?, "--hifi-accuracy")?;
                    config.hifi_error_rate = 1.0 - accuracy;
                }
                "--minimap-min-identity" => {
                    i += 1;
                    config.minimap_min_identity = parse_f64(
                        take_arg(&args, i, "minimap minimum identity")?,
                        "--minimap-min-identity",
                    )?;
                }
                "--minimap-min-align-len" => {
                    i += 1;
                    config.minimap_min_align_len = parse_usize(
                        take_arg(&args, i, "minimap minimum alignment length")?,
                        "--minimap-min-align-len",
                    )?;
                }
                "--paf-max-link-gap" => {
                    i += 1;
                    config.paf_max_link_gap = take_arg(&args, i, "PAF max link gap")?
                        .parse::<isize>()
                        .map_err(|_| "invalid --paf-max-link-gap".to_string())?;
                }
                "--skeleton-gfa" => {
                    i += 1;
                    config.skeleton_gfa = Some(PathBuf::from(take_arg(&args, i, "skeleton GFA")?));
                }
                "--skeleton-only" => {
                    config.skeleton_only = true;
                }
                "--skeleton-end-slop" => {
                    i += 1;
                    config.skeleton_end_slop = parse_usize(
                        take_arg(&args, i, "skeleton end slop")?,
                        "--skeleton-end-slop",
                    )?;
                }
                "--skeleton-min-link-support" => {
                    i += 1;
                    config.skeleton_min_link_support = parse_u32(
                        take_arg(&args, i, "skeleton minimum link support")?,
                        "--skeleton-min-link-support",
                    )?;
                }
                "--skeleton-min-link-ratio" => {
                    i += 1;
                    config.skeleton_min_link_ratio = parse_f64(
                        take_arg(&args, i, "skeleton minimum link ratio")?,
                        "--skeleton-min-link-ratio",
                    )?;
                }
                "--skeleton-rescue-gfa" => {
                    i += 1;
                    config.skeleton_rescue_gfa =
                        Some(PathBuf::from(take_arg(&args, i, "skeleton rescue GFA")?));
                }
                "--skeleton-rescue-link-support" => {
                    i += 1;
                    config.skeleton_rescue_link_support = parse_u32(
                        take_arg(&args, i, "skeleton rescue link support")?,
                        "--skeleton-rescue-link-support",
                    )?;
                }
                "--rounds" => {
                    i += 1;
                    config.rounds = parse_usize(take_arg(&args, i, "rounds")?, "--rounds")?;
                    overrides.rounds = true;
                }
                "-t" | "--threads" => {
                    i += 1;
                    config.threads = parse_usize(take_arg(&args, i, "threads")?, "--threads")?;
                }
                "--minimap2" => {
                    config.run_minimap2 = true;
                }
                "--help-advanced" => {
                    print_advanced_usage();
                    std::process::exit(0);
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

        apply_profiles(&mut config, &overrides);
        Ok(config)
    }
}

fn apply_profiles(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.rounds && config.rounds == 0 {
        config.rounds = match config.organelle {
            Some(OrganelleProfile::Mito) => 2,
            _ => 1,
        };
    }
    if config.organelle == Some(OrganelleProfile::Mito) && !overrides.numt_interference {
        config.numt_interference = NumtInterference::High;
    }

    match (config.organelle, config.numt_interference) {
        (Some(OrganelleProfile::Plastid), _) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 100;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 100;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.25;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 4;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 1000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.read_junction_links {
                config.read_junction_links = true;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (Some(OrganelleProfile::Mito), NumtInterference::Low) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 18;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 18;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.25;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 4;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 3000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (Some(OrganelleProfile::Mito), NumtInterference::High) | (None, NumtInterference::High) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 18;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 18;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.30;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 3;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 3000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (None, NumtInterference::Low) => {
            if config.rounds == 0 {
                config.rounds = 1;
            }
        }
    }
}

fn apply_common_organelle_defaults(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.k {
        config.k = 251;
    }
    if !overrides.s {
        config.s = 21;
    }
    if !overrides.dedup_kmer {
        config.dedup_kmer = 17;
    }
    if !overrides.containment_ratio {
        config.containment_ratio = 0.70;
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
         Common options:\n\
           --organelle plastid|mito     plastid defaults to 1 round; mito defaults to 2 rounds\n\
           --rounds INT                 override profile rounds\n\
           --numt-interference low|high mito profile strictness, mito default: high\n\
           -i, --pacbio-hifi FILE       input reads; may be repeated\n\
           -o DIR                       output directory\n\
           -t INT                       threads\n\
           --help-advanced              show tuning and skeleton-linking options\n\
         \n\
         Profile defaults:\n\
           plastid: k=251, s=21, coverage=100, clean 1-round graph\n\
           mito:    k=251, s=21, coverage=18, high NUMT caution, terminal skeleton + endpoint-link round\n\
         \n\
         Typical commands:\n\
           simple_draft_asm --organelle plastid -i data/plastid.fastq.gz -o result_plastid\n\
           simple_draft_asm --organelle mito -i data/mito.fastq.gz -o result_mito"
    );
}

fn print_advanced_usage() {
    eprintln!(
        "Advanced options supported by simple_draft_asm:\n\
         \n\
         Anchor selection:\n\
           --anchor syncmer|minimizer\n\
           --k INT\n\
           --s INT\n\
           --syncmer-pos INT\n\
           --window INT\n\
           --min-anchor-spacing INT\n\
         \n\
         First-round graph filtering:\n\
           --min-anchor-coverage INT\n\
           --min-edge-coverage INT\n\
           --min-branch-ratio FLOAT\n\
           --max-edges-per-state INT\n\
           --junction-rescue-support INT\n\
           --dedup-kmer INT\n\
           --containment-ratio FLOAT\n\
           --min-unitig-len INT\n\
           --min-tip-len INT\n\
           --min-link-support INT\n\
           --read-junction-links\n\
           --bidirectional-links\n\
         \n\
         HiFi and minimap2 support:\n\
           --hifi-error-rate FLOAT\n\
           --hifi-accuracy FLOAT\n\
           --minimap2\n\
           --minimap-min-identity FLOAT\n\
           --minimap-min-align-len INT\n\
           --paf-max-link-gap INT\n\
         \n\
         Skeleton second round:\n\
           --skeleton-gfa FILE\n\
           --skeleton-only\n\
           --skeleton-end-slop INT\n\
           --skeleton-min-link-support INT\n\
           --skeleton-min-link-ratio FLOAT\n\
           --skeleton-rescue-gfa FILE\n\
           --skeleton-rescue-link-support INT"
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

fn write_unitigs(config: &Config, unitigs: &[Unitig], filename: &str) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join(filename))?;
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
    filename: &str,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join(filename))?;
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

fn write_junctions(
    config: &Config,
    junctions: &HashMap<(usize, usize), u32>,
    filename: &str,
) -> io::Result<()> {
    let mut out_totals: HashMap<usize, u32> = HashMap::new();
    for (&(from, _to), &count) in junctions {
        *out_totals.entry(from).or_insert(0) += count;
    }

    let mut rows: Vec<_> = junctions.iter().collect();
    rows.sort_by_key(|(&(from, to), _)| (from, to));

    let mut out = File::create(config.out_dir.join(filename))?;
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
    full_graph: &FilteredGraph,
    unitigs: &[Unitig],
    full_unitigs: &[Unitig],
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
    writeln!(out, "full_edges\t{}", full_graph.edges.len())?;
    writeln!(out, "unitigs\t{}", unitigs.len())?;
    writeln!(
        out,
        "unitig_bases\t{}",
        unitigs
            .iter()
            .map(|unitig| unitig.sequence.len())
            .sum::<usize>()
    )?;
    writeln!(out, "full_unitigs\t{}", full_unitigs.len())?;
    writeln!(
        out,
        "full_unitig_bases\t{}",
        full_unitigs
            .iter()
            .map(|unitig| unitig.sequence.len())
            .sum::<usize>()
    )?;
    writeln!(
        out,
        "organelle\t{}",
        config
            .organelle
            .map(|o| o.as_str())
            .unwrap_or("unspecified")
    )?;
    writeln!(
        out,
        "numt_interference\t{}",
        config.numt_interference.as_str()
    )?;
    writeln!(out, "rounds\t{}", config.rounds)?;
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
    writeln!(out, "hifi_error_rate\t{:.5}", config.hifi_error_rate)?;
    writeln!(
        out,
        "expected_exact_anchor_survival\t{:.6}",
        expected_exact_survival(config.k, config.hifi_error_rate)
    )?;
    writeln!(
        out,
        "expected_exact_anchor_pair_survival\t{:.6}",
        expected_exact_survival(config.k, config.hifi_error_rate).powi(2)
    )?;
    writeln!(out, "minimap2\t{}", config.run_minimap2)?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    writeln!(out, "threads\t{}", config.threads)?;
    writeln!(out, "graph_gfa\tgraph.gfa")?;
    writeln!(out, "graph_full_gfa\tgraph.full.gfa")?;
    writeln!(
        out,
        "support_note\tRC:i is observed exact-anchor/read-walk support; long-k HiFi anchor dropout can underestimate true support, so use minimap2 depth/link TSVs when --minimap2 is enabled."
    )?;
    Ok(())
}

fn expected_exact_survival(k: usize, error_rate: f64) -> f64 {
    (1.0 - error_rate).clamp(0.0, 1.0).powi(k as i32)
}

#[derive(Debug, Clone)]
struct PafAln {
    qname: String,
    qstart: usize,
    qend: usize,
    strand: char,
    tname: String,
    tlen: usize,
    tstart: usize,
    tend: usize,
    matches: usize,
    alen: usize,
}

impl PafAln {
    fn identity(&self) -> f64 {
        if self.alen == 0 {
            0.0
        } else {
            self.matches as f64 / self.alen as f64
        }
    }
}

#[derive(Debug, Clone)]
struct SkeletonSegment {
    name: String,
    sequence: String,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct SkeletonLinkKey {
    from: String,
    from_orient: char,
    to: String,
    to_orient: char,
}

#[derive(Debug, Clone, Default)]
struct SkeletonLinkSupport {
    skeleton_support: u32,
    paf_support: u32,
    rescue_support: u32,
    out_ratio: f64,
    in_ratio: f64,
}

fn run_skeleton_linking(config: &Config) -> io::Result<()> {
    let skeleton_gfa = config
        .skeleton_gfa
        .as_ref()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "--skeleton-gfa is required"))?;
    let (segments, skeleton_links) = read_skeleton_gfa(skeleton_gfa)?;
    if segments.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} has no S records", skeleton_gfa.display()),
        ));
    }

    fs::create_dir_all(&config.out_dir)?;
    let skeleton_fasta = config.out_dir.join("skeleton.segments.fasta");
    write_skeleton_fasta(&skeleton_fasta, &segments)?;

    let paf_path = config.out_dir.join("reads_to_skeleton.paf");
    let log_path = config.out_dir.join("reads_to_skeleton.minimap2.log");
    run_minimap2_to_target(config, &skeleton_fasta, &paf_path, &log_path)?;

    let (depth, paf_links) = summarize_skeleton_paf(config, &segments, &paf_path)?;
    write_skeleton_depth(
        &config.out_dir.join("skeleton.depth.tsv"),
        &segments,
        &depth,
    )?;
    let mut merged_links = merge_skeleton_links(config, skeleton_links, paf_links);
    if let Some(rescue_gfa) = &config.skeleton_rescue_gfa {
        let rescue_links = read_gfa_links(rescue_gfa)?;
        add_component_rescue_links(config, &segments, &mut merged_links, rescue_links);
    }
    write_skeleton_links(&config.out_dir.join("skeleton.links.tsv"), &merged_links)?;
    write_skeleton_linked_gfa(
        config,
        &config.out_dir.join("skeleton.linked.gfa"),
        &segments,
        &depth,
        &merged_links,
    )?;
    write_skeleton_report(config, skeleton_gfa, &segments, &merged_links)?;
    Ok(())
}

fn read_skeleton_gfa(
    path: &Path,
) -> io::Result<(
    Vec<SkeletonSegment>,
    HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
)> {
    let reader = BufReader::new(File::open(path)?);
    let mut segments = Vec::new();
    let mut links: HashMap<SkeletonLinkKey, SkeletonLinkSupport> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        match fields[0] {
            "S" if fields.len() >= 3 => segments.push(SkeletonSegment {
                name: fields[1].to_string(),
                sequence: fields[2].to_string(),
            }),
            "L" if fields.len() >= 5 => {
                let mut support = 0u32;
                for field in &fields[5..] {
                    if let Some(value) = field.strip_prefix("RC:i:") {
                        support = value.parse().unwrap_or(0);
                    }
                }
                let key = SkeletonLinkKey {
                    from: fields[1].to_string(),
                    from_orient: fields[2].chars().next().unwrap_or('+'),
                    to: fields[3].to_string(),
                    to_orient: fields[4].chars().next().unwrap_or('+'),
                };
                links.entry(key).or_default().skeleton_support = support;
            }
            _ => {}
        }
    }
    Ok((segments, links))
}

fn read_gfa_links(path: &Path) -> io::Result<HashMap<SkeletonLinkKey, u32>> {
    let reader = BufReader::new(File::open(path)?);
    let mut links: HashMap<SkeletonLinkKey, u32> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 || fields[0] != "L" {
            continue;
        }
        let mut support = 0u32;
        for field in &fields[5..] {
            if let Some(value) = field.strip_prefix("RC:i:") {
                support = value.parse().unwrap_or(0);
            }
        }
        let key = SkeletonLinkKey {
            from: fields[1].to_string(),
            from_orient: fields[2].chars().next().unwrap_or('+'),
            to: fields[3].to_string(),
            to_orient: fields[4].chars().next().unwrap_or('+'),
        };
        links
            .entry(key)
            .and_modify(|old| *old = (*old).max(support))
            .or_insert(support);
    }
    Ok(links)
}

fn add_component_rescue_links(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    rescue_links: HashMap<SkeletonLinkKey, u32>,
) {
    let components = skeleton_components(config, segments, links);
    for (key, support) in rescue_links {
        if support < config.skeleton_rescue_link_support {
            continue;
        }
        let Some(&from_component) = components.get(&key.from) else {
            continue;
        };
        let Some(&to_component) = components.get(&key.to) else {
            continue;
        };
        if from_component == to_component {
            continue;
        }
        links.entry(key).or_default().rescue_support = support;
    }
}

fn skeleton_components(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, usize> {
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();
    for segment in segments {
        adj.entry(segment.name.clone()).or_default();
    }
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        adj.entry(key.from.clone())
            .or_default()
            .push(key.to.clone());
        adj.entry(key.to.clone())
            .or_default()
            .push(key.from.clone());
    }
    let mut component = HashMap::new();
    let mut next_component = 0usize;
    for segment in segments {
        if component.contains_key(&segment.name) {
            continue;
        }
        next_component += 1;
        let mut stack = vec![segment.name.clone()];
        component.insert(segment.name.clone(), next_component);
        while let Some(node) = stack.pop() {
            if let Some(neighbors) = adj.get(&node) {
                for neighbor in neighbors {
                    if !component.contains_key(neighbor) {
                        component.insert(neighbor.clone(), next_component);
                        stack.push(neighbor.clone());
                    }
                }
            }
        }
    }
    component
}

fn write_skeleton_fasta(path: &Path, segments: &[SkeletonSegment]) -> io::Result<()> {
    let mut out = File::create(path)?;
    for segment in segments {
        writeln!(out, ">{}", segment.name)?;
        write_wrapped(&mut out, &segment.sequence, 80)?;
    }
    Ok(())
}

fn run_minimap2_to_target(
    config: &Config,
    target: &Path,
    paf_path: &Path,
    log_path: &Path,
) -> io::Result<()> {
    let paf = File::create(paf_path)?;
    let mut command = Command::new("minimap2");
    command
        .arg("-x")
        .arg("map-hifi")
        .arg("-t")
        .arg(config.threads.to_string())
        .arg(target);
    for read in &config.reads {
        command.arg(read);
    }
    let output = command
        .stdout(Stdio::from(paf))
        .stderr(Stdio::piped())
        .output();
    match output {
        Ok(output) if output.status.success() => {
            fs::write(log_path, &output.stderr)?;
            Ok(())
        }
        Ok(output) => {
            fs::write(log_path, &output.stderr)?;
            Err(io::Error::other(format!(
                "minimap2 failed; see {}",
                log_path.display()
            )))
        }
        Err(err) if err.kind() == io::ErrorKind::NotFound => Err(io::Error::new(
            io::ErrorKind::NotFound,
            "minimap2 not found in PATH",
        )),
        Err(err) => Err(err),
    }
}

fn summarize_skeleton_paf(
    config: &Config,
    segments: &[SkeletonSegment],
    paf_path: &Path,
) -> io::Result<(
    HashMap<String, (usize, usize)>,
    HashMap<SkeletonLinkKey, u32>,
)> {
    let lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let mut depth: HashMap<String, (usize, usize)> = HashMap::new();
    let mut by_read: HashMap<String, Vec<PafAln>> = HashMap::new();
    let reader = BufReader::new(File::open(paf_path)?);
    for line in reader.lines() {
        let line = line?;
        let Some(aln) = parse_paf_line(&line) else {
            continue;
        };
        if aln.identity() < config.minimap_min_identity || aln.alen < config.minimap_min_align_len {
            continue;
        }
        if !lengths.contains_key(&aln.tname) {
            continue;
        }
        let entry = depth.entry(aln.tname.clone()).or_insert((0, 0));
        entry.0 += aln.tend.saturating_sub(aln.tstart);
        entry.1 += 1;
        by_read.entry(aln.qname.clone()).or_default().push(aln);
    }

    let mut links: HashMap<SkeletonLinkKey, u32> = HashMap::new();
    for (_read, mut alns) in by_read {
        alns.sort_by(|a, b| {
            a.qstart
                .cmp(&b.qstart)
                .then_with(|| b.alen.cmp(&a.alen))
                .then_with(|| {
                    b.identity()
                        .partial_cmp(&a.identity())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        let chain = paf_non_overlapping_chain(alns, config.paf_max_link_gap);
        for pair in chain.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            let gap = right.qstart as isize - left.qend as isize;
            if gap.abs() > config.paf_max_link_gap {
                continue;
            }
            let Some(from_orient) = paf_exit_orient(left, config.skeleton_end_slop) else {
                continue;
            };
            let Some(to_orient) = paf_entry_orient(right, config.skeleton_end_slop) else {
                continue;
            };
            if left.tname == right.tname && from_orient == to_orient {
                continue;
            }
            let key = SkeletonLinkKey {
                from: left.tname.clone(),
                from_orient,
                to: right.tname.clone(),
                to_orient,
            };
            *links.entry(key).or_insert(0) += 1;
        }
    }
    Ok((depth, links))
}

fn paf_non_overlapping_chain(mut alns: Vec<PafAln>, max_gap: isize) -> Vec<PafAln> {
    let mut chain: Vec<PafAln> = Vec::new();
    alns.sort_by(|a, b| {
        a.qstart
            .cmp(&b.qstart)
            .then_with(|| b.alen.cmp(&a.alen))
            .then_with(|| {
                b.identity()
                    .partial_cmp(&a.identity())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    for aln in alns {
        if let Some(last) = chain.last() {
            let overlap = last.qend as isize - aln.qstart as isize;
            if overlap > max_gap {
                continue;
            }
        }
        chain.push(aln);
    }
    chain
}

fn paf_exit_orient(aln: &PafAln, slop: usize) -> Option<char> {
    match aln.strand {
        '+' if aln.tlen.saturating_sub(aln.tend) <= slop => Some('+'),
        '-' if aln.tstart <= slop => Some('-'),
        _ => None,
    }
}

fn paf_entry_orient(aln: &PafAln, slop: usize) -> Option<char> {
    match aln.strand {
        '+' if aln.tstart <= slop => Some('+'),
        '-' if aln.tlen.saturating_sub(aln.tend) <= slop => Some('-'),
        _ => None,
    }
}

fn merge_skeleton_links(
    config: &Config,
    mut links: HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    paf_links: HashMap<SkeletonLinkKey, u32>,
) -> HashMap<SkeletonLinkKey, SkeletonLinkSupport> {
    for (key, support) in paf_links {
        links.entry(key).or_default().paf_support += support;
    }

    let mut out_max: HashMap<(String, char), u32> = HashMap::new();
    let mut in_max: HashMap<(String, char), u32> = HashMap::new();
    for (key, support) in &links {
        if support.paf_support == 0 {
            continue;
        }
        out_max
            .entry((key.from.clone(), key.from_orient))
            .and_modify(|max| *max = (*max).max(support.paf_support))
            .or_insert(support.paf_support);
        in_max
            .entry((key.to.clone(), key.to_orient))
            .and_modify(|max| *max = (*max).max(support.paf_support))
            .or_insert(support.paf_support);
    }
    for (key, support) in &mut links {
        let out_best = out_max
            .get(&(key.from.clone(), key.from_orient))
            .copied()
            .unwrap_or(support.paf_support.max(1));
        let in_best = in_max
            .get(&(key.to.clone(), key.to_orient))
            .copied()
            .unwrap_or(support.paf_support.max(1));
        support.out_ratio = support.paf_support as f64 / out_best.max(1) as f64;
        support.in_ratio = support.paf_support as f64 / in_best.max(1) as f64;
    }

    if config.bidirectional_links {
        let mut extra = Vec::new();
        for (key, support) in &links {
            let rc = SkeletonLinkKey {
                from: key.to.clone(),
                from_orient: flip_orient(key.to_orient),
                to: key.from.clone(),
                to_orient: flip_orient(key.from_orient),
            };
            if !links.contains_key(&rc) {
                extra.push((rc, support.clone()));
            }
        }
        for (key, support) in extra {
            links.insert(key, support);
        }
    }
    links
}

fn write_skeleton_depth(
    path: &Path,
    segments: &[SkeletonSegment],
    depth: &HashMap<String, (usize, usize)>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "segment\tlength\taligned_bases\tmean_depth\talignments"
    )?;
    for segment in segments {
        let (bases, count) = depth.get(&segment.name).copied().unwrap_or((0, 0));
        let mean = if segment.sequence.is_empty() {
            0.0
        } else {
            bases as f64 / segment.sequence.len() as f64
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{:.3}\t{}",
            segment.name,
            segment.sequence.len(),
            bases,
            mean,
            count
        )?;
    }
    Ok(())
}

fn write_skeleton_links(
    path: &Path,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "from\tfrom_orient\tto\tto_orient\tskeleton_support\tpaf_support\trescue_support\tout_ratio\tin_ratio"
    )?;
    let mut rows: Vec<_> = links.iter().collect();
    rows.sort_by(|a, b| {
        b.1.paf_support
            .cmp(&a.1.paf_support)
            .then_with(|| b.1.skeleton_support.cmp(&a.1.skeleton_support))
            .then_with(|| b.1.rescue_support.cmp(&a.1.rescue_support))
            .then_with(|| a.0.cmp(b.0))
    });
    for (key, support) in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            support.skeleton_support,
            support.paf_support,
            support.rescue_support,
            support.out_ratio,
            support.in_ratio
        )?;
    }
    Ok(())
}

fn write_skeleton_linked_gfa(
    config: &Config,
    path: &Path,
    segments: &[SkeletonSegment],
    depth: &HashMap<String, (usize, usize)>,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "H\tVN:Z:1.0\tPG:Z:simple_draft_asm\tST:Z:skeleton_linked"
    )?;
    for segment in segments {
        let (aligned_bases, alignments) = depth.get(&segment.name).copied().unwrap_or((0, 0));
        let mean_depth = if segment.sequence.is_empty() {
            0.0
        } else {
            aligned_bases as f64 / segment.sequence.len() as f64
        };
        writeln!(
            out,
            "S\t{}\t{}\tLN:i:{}\tDP:f:{:.3}\tAB:i:{}\tAC:i:{}",
            segment.name,
            segment.sequence,
            segment.sequence.len(),
            mean_depth,
            aligned_bases,
            alignments
        )?;
    }
    let mut rows: Vec<_> = links.iter().collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for (key, support) in rows {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        let ratio = support.out_ratio.min(support.in_ratio);
        let rc = support
            .skeleton_support
            .max(support.paf_support)
            .max(support.rescue_support);
        let source = match (support.skeleton_support > 0, support.paf_support > 0) {
            (true, true) if support.rescue_support > 0 => "skeleton+paf+rescue",
            (true, true) => "skeleton+paf",
            (true, false) if support.rescue_support > 0 => "skeleton+rescue",
            (true, false) => "skeleton",
            (false, true) if support.rescue_support > 0 => "paf+rescue",
            (false, true) => "paf",
            (false, false) if support.rescue_support > 0 => "rescue",
            (false, false) => "none",
        };
        writeln!(
            out,
            "L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\tSK:i:{}\tPA:i:{}\tRS:i:{}\tLR:f:{:.4}\tSC:Z:{}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            rc,
            support.skeleton_support,
            support.paf_support,
            support.rescue_support,
            ratio,
            source
        )?;
    }
    Ok(())
}

fn write_skeleton_report(
    config: &Config,
    skeleton_gfa: &Path,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("skeleton.report.txt"))?;
    let kept = links
        .values()
        .filter(|support| keep_skeleton_link(config, support))
        .count();
    writeln!(out, "skeleton_gfa\t{}", skeleton_gfa.display())?;
    writeln!(out, "segments\t{}", segments.len())?;
    writeln!(out, "links_total\t{}", links.len())?;
    writeln!(out, "links_kept\t{}", kept)?;
    writeln!(out, "skeleton_end_slop\t{}", config.skeleton_end_slop)?;
    writeln!(
        out,
        "skeleton_min_link_support\t{}",
        config.skeleton_min_link_support
    )?;
    writeln!(
        out,
        "skeleton_min_link_ratio\t{}",
        config.skeleton_min_link_ratio
    )?;
    if let Some(rescue_gfa) = &config.skeleton_rescue_gfa {
        writeln!(out, "skeleton_rescue_gfa\t{}", rescue_gfa.display())?;
    }
    writeln!(
        out,
        "skeleton_rescue_link_support\t{}",
        config.skeleton_rescue_link_support
    )?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    writeln!(out, "output_gfa\tskeleton.linked.gfa")?;
    Ok(())
}

fn keep_skeleton_link(config: &Config, support: &SkeletonLinkSupport) -> bool {
    support.skeleton_support >= config.min_link_support
        || (support.paf_support >= config.skeleton_min_link_support
            && support.out_ratio.min(support.in_ratio) >= config.skeleton_min_link_ratio)
        || support.rescue_support >= config.skeleton_rescue_link_support
}

fn flip_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn run_minimap2(config: &Config) -> io::Result<()> {
    let jobs = [
        (
            "dominant",
            config.out_dir.join("unitigs.fasta"),
            config.out_dir.join("read_to_unitigs.paf"),
            config.out_dir.join("read_to_unitigs.minimap2.log"),
            config.out_dir.join("depth.minimap2.tsv"),
            config.out_dir.join("junctions.minimap2.tsv"),
        ),
        (
            "full",
            config.out_dir.join("unitigs.full.fasta"),
            config.out_dir.join("read_to_unitigs.full.paf"),
            config.out_dir.join("read_to_unitigs.full.minimap2.log"),
            config.out_dir.join("depth.minimap2.full.tsv"),
            config.out_dir.join("junctions.minimap2.full.tsv"),
        ),
    ];

    for (label, unitigs, paf_path, log_path, depth_path, junction_path) in jobs {
        if fasta_is_empty(&unitigs)? {
            continue;
        }
        let paf = File::create(&paf_path)?;
        let mut command = Command::new("minimap2");
        command
            .arg("-x")
            .arg("map-hifi")
            .arg("-t")
            .arg(config.threads.to_string())
            .arg(&unitigs);
        for read in &config.reads {
            command.arg(read);
        }
        let output = command
            .stdout(Stdio::from(paf))
            .stderr(Stdio::piped())
            .output();
        match output {
            Ok(output) if output.status.success() => {
                fs::write(&log_path, &output.stderr)?;
                summarize_paf(config, &paf_path, &unitigs, &depth_path, &junction_path)?;
            }
            Ok(output) => {
                fs::write(&log_path, &output.stderr)?;
                return Err(io::Error::other(format!(
                    "minimap2 failed for {label}; see {}",
                    log_path.display()
                )));
            }
            Err(err) if err.kind() == io::ErrorKind::NotFound => {
                let message = "warning: --minimap2 requested, but minimap2 was not found in PATH\n";
                fs::write(&log_path, message)?;
                eprint!("{message}");
            }
            Err(err) => return Err(err),
        }
    }
    Ok(())
}

fn fasta_is_empty(path: &Path) -> io::Result<bool> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if line?.starts_with('>') {
            return Ok(false);
        }
    }
    Ok(true)
}

fn summarize_paf(
    config: &Config,
    paf_path: &Path,
    target_fasta: &Path,
    depth_path: &Path,
    junction_path: &Path,
) -> io::Result<()> {
    let mut target_lengths = read_fasta_lengths(target_fasta)?;
    let mut depth_bases: HashMap<String, usize> = HashMap::new();
    let mut alignment_counts: HashMap<String, usize> = HashMap::new();
    let mut by_read: HashMap<String, Vec<PafAln>> = HashMap::new();

    let reader = BufReader::new(File::open(paf_path)?);
    for line in reader.lines() {
        let line = line?;
        let Some(aln) = parse_paf_line(&line) else {
            continue;
        };
        if aln.identity() < config.minimap_min_identity || aln.alen < config.minimap_min_align_len {
            continue;
        }
        target_lengths.entry(aln.tname.clone()).or_insert(aln.tlen);
        *depth_bases.entry(aln.tname.clone()).or_insert(0) += aln.tend.saturating_sub(aln.tstart);
        *alignment_counts.entry(aln.tname.clone()).or_insert(0) += 1;
        by_read.entry(aln.qname.clone()).or_default().push(aln);
    }

    let mut depth_out = File::create(depth_path)?;
    writeln!(
        depth_out,
        "unitig\tlength\taligned_bases\tmean_depth\talignments"
    )?;
    let mut names: Vec<_> = target_lengths.keys().cloned().collect();
    names.sort();
    for name in names {
        let len = target_lengths.get(&name).copied().unwrap_or(0);
        let bases = depth_bases.get(&name).copied().unwrap_or(0);
        let count = alignment_counts.get(&name).copied().unwrap_or(0);
        let depth = if len == 0 {
            0.0
        } else {
            bases as f64 / len as f64
        };
        writeln!(depth_out, "{name}\t{len}\t{bases}\t{depth:.3}\t{count}")?;
    }

    let mut links: HashMap<(String, char, String, char), usize> = HashMap::new();
    for (_read, mut alns) in by_read {
        alns.sort_by(|a, b| {
            a.qstart
                .cmp(&b.qstart)
                .then_with(|| b.alen.cmp(&a.alen))
                .then_with(|| {
                    b.identity()
                        .partial_cmp(&a.identity())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        let mut chain: Vec<PafAln> = Vec::new();
        for aln in alns {
            if let Some(last) = chain.last() {
                let overlap = last.qend as isize - aln.qstart as isize;
                if overlap > config.paf_max_link_gap {
                    continue;
                }
            }
            chain.push(aln);
        }
        for pair in chain.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            if left.tname == right.tname && left.strand == right.strand {
                continue;
            }
            let gap = right.qstart as isize - left.qend as isize;
            if gap.abs() <= config.paf_max_link_gap {
                *links
                    .entry((
                        left.tname.clone(),
                        left.strand,
                        right.tname.clone(),
                        right.strand,
                    ))
                    .or_insert(0) += 1;
            }
        }
    }

    let mut junction_out = File::create(junction_path)?;
    writeln!(
        junction_out,
        "from_unitig\tfrom_orient\tto_unitig\tto_orient\tcount"
    )?;
    let mut rows: Vec<_> = links.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    for ((from, from_orient, to, to_orient), count) in rows {
        writeln!(
            junction_out,
            "{from}\t{from_orient}\t{to}\t{to_orient}\t{count}"
        )?;
    }
    Ok(())
}

fn parse_paf_line(line: &str) -> Option<PafAln> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return None;
    }
    Some(PafAln {
        qname: fields[0].to_string(),
        qstart: fields[2].parse().ok()?,
        qend: fields[3].parse().ok()?,
        strand: fields[4].chars().next()?,
        tname: fields[5].to_string(),
        tlen: fields[6].parse().ok()?,
        tstart: fields[7].parse().ok()?,
        tend: fields[8].parse().ok()?,
        matches: fields[9].parse().ok()?,
        alen: fields[10].parse().ok()?,
    })
}

fn read_fasta_lengths(path: &Path) -> io::Result<HashMap<String, usize>> {
    let reader = BufReader::new(File::open(path)?);
    let mut lengths = HashMap::new();
    let mut current = None;
    let mut len = 0usize;
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current.take() {
                lengths.insert(name, len);
            }
            current = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string(),
            );
            len = 0;
        } else {
            len += line.trim().len();
        }
    }
    if let Some(name) = current {
        lengths.insert(name, len);
    }
    Ok(lengths)
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
