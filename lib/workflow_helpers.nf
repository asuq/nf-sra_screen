// Shared workflow helper functions for nf-sra_screen entrypoints and local subworkflows.

/*
 * Return all implemented binners in their stable execution order.
 */
def allImplementedBinners() {
  return ['metabat', 'semibin', 'rosella', 'comebin', 'vamb', 'lorbin']
}


/*
 * Return the binners compatible with the supplied read type.
 */
def compatibleBinnersForReadType(readType) {
  def binners = ['metabat', 'semibin', 'rosella', 'comebin', 'vamb']
  if (readType?.toString()?.trim()?.equalsIgnoreCase('hifi')) {
    binners = binners + ['lorbin']
  }
  return binners
}


/*
 * Parse and validate a comma-separated binning tool selection.
 */
def parseBinnerSelection(raw_selection_value) {
  def normalised_selection_text = raw_selection_value == null ? 'all-compatible' : raw_selection_value.toString().trim().toLowerCase()
  if (!normalised_selection_text) {
    error "--binners must include at least one tool"
  }

  if (normalised_selection_text in ['auto', 'all-compatible']) {
    return [mode: 'auto', tools: allImplementedBinners()]
  }

  return [
    mode: 'explicit',
    tools: parseToolSelection(
      normalised_selection_text,
      'all-compatible',
      allImplementedBinners() as Set,
      [] as Set,
      'binners'
    )
  ]
}


/*
 * Return the effective binners for one sample, or fail on invalid LorBin use.
 */
def binnersForSample(selection, sra, srr, readType) {
  def readTypeLc = readType?.toString()?.trim()?.toLowerCase()
  def mode = selection.mode ?: selection.binnerMode
  def tools = selection.tools ?: selection.binners
  def selected_binners = mode == 'auto'
    ? compatibleBinnersForReadType(readTypeLc)
    : tools

  if ('lorbin' in selected_binners && readTypeLc != 'hifi') {
    error "LorBin only supports hifi reads; sample ${sra}:${srr} has read_type ${readType}. Remove lorbin from --binners for this run."
  }

  return selected_binners
}


/*
 * Test whether a comma-separated binner list includes one tool.
 */
def binnerCsvContains(binner_csv, tool) {
  return binner_csv.toString().split(',').collect { token -> token.trim() }.contains(tool)
}


/*
 * Count tools in a comma-separated binner list.
 */
def binnerCsvSize(binner_csv) {
  return binner_csv.toString().split(',').collect { token -> token.trim() }.findAll { token -> token }.size()
}


/*
 * Parse and validate a comma-separated tool selection.
 */
def parseToolSelection(raw_selection_value, defaultValue, allowed_tools, planned_tools, parameter_name) {
  def normalised_selection_text = raw_selection_value == null ? defaultValue : raw_selection_value.toString()
  def selected_tools = normalised_selection_text
    .split(',')
    .collect { token -> token.trim().toLowerCase() }
    .findAll { token -> token }
    .unique()

  if (!selected_tools) {
    error "--${parameter_name} must include at least one tool"
  }

  def planned_but_unimplemented_tools = selected_tools.findAll { tool_name -> tool_name in planned_tools }
  if (planned_but_unimplemented_tools) {
    error "--${parameter_name} includes planned tool(s) not implemented yet: ${planned_but_unimplemented_tools.join(', ')}"
  }

  def unsupported_tools = selected_tools.findAll { tool_name -> !(tool_name in allowed_tools) }
  if (unsupported_tools) {
    error "--${parameter_name} includes unsupported tool(s): ${unsupported_tools.join(', ')}"
  }

  selected_tools
}


/*
 * Validate the selected SemiBin2 pretrained environment.
 */
def validateSemibinEnvironment(raw_selection_value) {
  def selected_environment = (raw_selection_value ?: 'global').toString().trim().toLowerCase()
  def allowed_environments = [
    'human_gut',
    'dog_gut',
    'ocean',
    'soil',
    'cat_gut',
    'human_oral',
    'mouse_gut',
    'pig_gut',
    'built_environment',
    'wastewater',
    'chicken_caecum',
    'global'
  ] as Set

  if (!selected_environment) {
    error "--semibin_environment must not be empty"
  }

  if (!(selected_environment in allowed_environments)) {
    error "--semibin_environment includes unsupported environment: ${selected_environment}"
  }

  selected_environment
}


/*
 * Validate reserved binning syntax before workflow construction.
 */
def validateBinningOptions() {
  def use_gpu_binners = params.gpu?.toString()?.toBoolean() ?: false
  def planned_tools = [] as Set
  def binnerSelection = parseBinnerSelection(params.binners)
  def refiners = parseToolSelection(
    params.refiners,
    'binette',
    ['dastool', 'binette'] as Set,
    planned_tools,
    'refiners'
  )

  if ('binette' in refiners && !params.checkm2_db) {
    error "--checkm2_db is required when --refiners includes binette"
  }

  [
    binners: binnerSelection.tools,
    binnerMode: binnerSelection.mode,
    refiners: refiners,
    gpu: use_gpu_binners,
    semibinEnvironment: validateSemibinEnvironment(params.semibin_environment)
  ]
}


/*
 * Return the canonical assembler name for a user-provided tool token.
 */
def canonicalAssemblerName(name) {
  def value = (name ?: '').toString().trim().toLowerCase()
  def aliases = [
    spades: 'metaspades',
    flye: 'metaflye'
  ]
  aliases.get(value, value)
}


/*
 * Parse and validate the assembler selection CLI parameter.
 */
def selectedAssemblerTokens() {
  def raw_selection_value = params.assembler ?: params.assemblers ?: 'auto'
  def rawTokens = raw_selection_value.toString()
    .split(/\s*,\s*/)
    .collect { token -> canonicalAssemblerName(token) }
    .findAll { token -> token }
    .unique()

  def tokens = rawTokens ?: ['auto']
  def supported = ['metaspades', 'unicycler', 'metaflye', 'myloasm'] as Set
  def keywords = ['auto', 'all'] as Set
  def unsupported = tokens.findAll { token -> !(token in supported) && !(token in keywords) }

  if (unsupported) {
    error "Unsupported assembler(s): ${unsupported.join(', ')}"
  }
  if ((tokens.contains('auto') || tokens.contains('all')) && tokens.size() > 1) {
    error "--assemblers cannot combine auto/all with explicit assembler names"
  }

  tokens
}


/*
 * Return all supported assemblers for a read type.
 */
def supportedAssemblersForReadType(readType) {
  def supportedByReadType = [
    short:    ['metaspades', 'unicycler'],
    nanopore: ['metaflye'],
    pacbio:   ['metaflye'],
    hifi:     ['metaflye', 'myloasm']
  ]
  supportedByReadType[(readType ?: '').toString().toLowerCase()] ?: []
}


/*
 * Return the default assembler for a read type when --assemblers auto is used.
 */
def defaultAssemblerForReadType(readType) {
  def defaults = [
    short:    'metaspades',
    nanopore: 'metaflye',
    pacbio:   'metaflye',
    hifi:     'myloasm'
  ]
  defaults[(readType ?: '').toString().toLowerCase()]
}


/*
 * Resolve the selected compatible assemblers for one read type.
 */
def assemblersForReadType(readType) {
  def tokens = selectedAssemblerTokens()
  def supported = supportedAssemblersForReadType(readType)

  if (tokens == ['auto']) {
    def defaultAssembler = defaultAssemblerForReadType(readType)
    return defaultAssembler ? [defaultAssembler] : []
  }
  if (tokens == ['all']) {
    return supported
  }

  tokens.findAll { token -> token in supported }
}


/*
 * Return true when selected outputs need per-assembler publish directories.
 */
def useAssemblerSubdirectories() {
  def tokens = selectedAssemblerTokens()
  tokens == ['all'] || tokens.size() > 1
}


/*
 * Return the publish directory for assembler-specific outputs.
 */
def assemblerPublishDir(sra, srr, assembler) {
  def base = "${params.outdir}/${sra}/${srr}"
  useAssemblerSubdirectories() ? "${base}/${assembler}" : base
}


/*
 * Parse and validate one binning samplesheet row.
 */
def normaliseBinningRow(row) {
  def sampleRaw = (row['sample'] ?: '').trim()
  def srrRaw = (row['srr'] ?: '').trim()
  def readTypeRaw = (row['read_type'] ?: '').trim().toLowerCase()
  def readsRaw = (row['reads'] ?: '').trim()
  def assemblyRaw = (row['assembly_fasta'] ?: '').trim()

  if (!sampleRaw || !assemblyRaw) {
    error "Invalid binning.tsv row with missing sample or assembly_fasta: ${row}"
  }

  def hasReads = !!readsRaw
  def hasSrr = !!srrRaw
  if (hasReads == hasSrr) {
    error "Sample '${sampleRaw}' must provide exactly one of reads or srr"
  }

  def assemblyFasta = file(assemblyRaw, checkIfExists: true)

  if (hasSrr) {
    return [
      mode: 'srr',
      sra: sampleRaw,
      srr: srrRaw,
      assembly_fasta: assemblyFasta
    ]
  }

  if (!readTypeRaw) {
    error "Sample '${sampleRaw}' is missing read_type for local reads"
  }

  def supportedReadTypes = ['short', 'nanopore', 'pacbio', 'hifi'] as Set
  if (!(readTypeRaw in supportedReadTypes)) {
    error "Unsupported read_type '${readTypeRaw}' for sample '${sampleRaw}'"
  }

  def readFiles = readsRaw
    .split(/\s*,\s*/)
    .findAll { token -> token }
    .collect { read_path -> file(read_path, checkIfExists: true) }

  if (!readFiles) {
    error "No read files were provided for sample '${sampleRaw}'"
  }

  if (readTypeRaw == 'short' && readFiles.size() > 2) {
    error "Sample '${sampleRaw}' has ${readFiles.size()} short-read FASTQs; expected one or two"
  }

  return [
    mode: 'local',
    sra: sampleRaw,
    srr: sampleRaw,
    platform: 'UNKNOWN',
    model: 'UNKNOWN',
    strategy: 'UNKNOWN',
    read_type: readTypeRaw,
    assembler: 'provided',
    assembly_fasta: assemblyFasta,
    reads: readFiles
  ]
}


/*
 * Parse the resolved SRR metadata emitted by the resolver helper.
 */
def parseResolvedMetadata(sra, srr, resolvedTsv) {
  def fields = file(resolvedTsv).text
    .readLines()
    .find { value -> value?.trim() }
    ?.split(/\t/, -1)
    ?.collect { token -> token.trim() }

  if (!fields || fields.size() != 4) {
    error "Invalid resolved metadata for sample '${sra}' run '${srr}'"
  }

  def (platform, model, strategy, readType) = fields
  if (!platform || !model || !strategy || !readType) {
    error "Incomplete resolved metadata for sample '${sra}' run '${srr}'"
  }

  tuple(sra, srr, platform, model, strategy, readType)
}
