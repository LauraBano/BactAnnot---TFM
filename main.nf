nextflow.enable.dsl = 2

/*
 * Flujo general de BactAnnot
 *
 * FASTA -> Bakta -> proteínas FAA
 *              -> eggNOG / InterProScan / KofamScan / PathoFact
 *              -> normalización -> integración -> consenso y concordancia
 */

// Carga los procesos que se van a conectar dentro del workflow
include { BAKTA } from './modules/local/bakta'
include { PARSE_BAKTA } from './modules/local/parse_bakta'
include { EGGNOG } from './modules/local/eggnog'
include { PARSE_EGGNOG } from './modules/local/parse_eggnog'
include { INTERPROSCAN } from './modules/local/interproscan'
include { PARSE_INTERPROSCAN } from './modules/local/parse_interproscan'
include { KOFAMSCAN } from './modules/local/kofamscan'
include { PARSE_KOFAMSCAN } from './modules/local/parse_kofamscan'
include { PATHOFACT } from './modules/local/pathofact'
include { PARSE_PATHOFACT } from './modules/local/parse_pathofact'
include { INTEGRATE_ANNOTATIONS } from './modules/local/integrate_annotations'
include { INTEGRATE_EVIDENCE } from './modules/local/integrate_evidence'
include { BUILD_CONSENSUS } from './modules/local/build_consensus'
include { BUILD_CONCORDANCE } from './modules/local/build_concordance'
include { BUILD_REPORT } from './modules/local/build_report'


// BactAnnot acepta un samplesheet o un único FASTA y normaliza ambos al mismo formato interno
params.input     = null
params.fasta     = null
params.sample_id = null
params.skip_report = false


// Rutas externas necesarias para ejecutar las herramientas de anotación
params.bakta_db  = null
params.eggnog_db = null
params.interproscan_path = null
params.kofam_db = null
params.pathofact_path = null
params.pathofact_db = null
params.pathofact_conda_prefix = null


workflow {

    main:

    // Comprueba que la entrada y las bases de datos necesarias estén bien configuradas
    if (
        !params.input &&
        !params.fasta &&
        !params.sample_id
    ) {
        error """
        No input was provided.

        Use either a CSV samplesheet:

            nextflow run main.nf \\
                --input samplesheet.csv \\
                -profile local

        or a single FASTA file:

            nextflow run main.nf \\
                --fasta genome.fna \\
                --sample_id sample_name \\
                -profile local
        """.stripIndent()
    }


    if (
        params.input &&
        (params.fasta || params.sample_id)
    ) {
        error """
        Ambiguous input.

        Provide either:

            --input samplesheet.csv

        or:

            --fasta genome.fna \\
            --sample_id sample_name

        Do not provide both input modes simultaneously.
        """.stripIndent()
    }


    if (
        !params.input &&
        (!params.fasta || !params.sample_id)
    ) {
        error """
        Single-sample mode requires both:

            --fasta genome.fna
            --sample_id sample_name
        """.stripIndent()
    }


    if (!params.bakta_db) {
        error """
        Bakta database path is not configured.

        Configure params.bakta_db in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.eggnog_db) {
        error """
        eggNOG-mapper database path is not configured.

        Configure params.eggnog_db in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.kofam_db) {
        error """
        KOfam database path is not configured.

        Configure params.kofam_db in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.interproscan_path) {
        error """
        InterProScan executable path is not configured.

        Configure params.interproscan_path in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.pathofact_path) {
        error """
        PathoFact2 installation path is not configured.

        Configure params.pathofact_path in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.pathofact_db) {
        error """
        PathoFact2 database path is not configured.

        Configure params.pathofact_db in:

            conf/local.config
        """.stripIndent()
    }


    if (!params.pathofact_conda_prefix) {
        error """
        PathoFact2 Snakemake Conda prefix is not configured.

        Configure params.pathofact_conda_prefix in:

            conf/local.config
        """.stripIndent()
    }


    // Si usamos un samplesheet, lee cada fila y la convierte en tuple(sample_id, fasta)
    if (params.input) {

        samples_ch = channel
            .fromPath(
                params.input,
                checkIfExists: true
            )
            .splitCsv(
                header: true,
                strip: true
            )
            .map { row ->

                // Comprueba que el CSV tenga las columnas necesarias
                if (
                    !row.containsKey('sample_id') ||
                    !row.containsKey('fasta')
                ) {
                    error """
                    Invalid samplesheet header.

                    Required columns:

                        sample_id,fasta
                    """.stripIndent()
                }

                // Comprueba que cada muestra tenga un identificador válido
                def sample_id = row.sample_id?.trim()
                def fasta_value = row.fasta?.trim()


                // Comprueba que cada muestra tenga un identificador válido
                if (!sample_id) {
                    error "Samplesheet contains an empty sample_id."
                }


                if (
                    !sample_id.matches(
                        '[A-Za-z0-9][A-Za-z0-9_.-]*'
                    )
                ) {
                    error """
                    Invalid sample_id: ${sample_id}

                    Allowed characters:

                        A-Z a-z 0-9 _ - .

                    The sample_id must start with a letter or number.
                    """.stripIndent()
                }


                // Comprueba la ruta del archivo FASTA
                if (!fasta_value) {
                    error """
                    Missing FASTA path for sample:

                        ${sample_id}
                    """.stripIndent()
                }

                // Convierte la ruta del FASTA en un archivo que Nextflow pueda manejar
                def fasta = file(
                    fasta_value,
                    checkIfExists: true
                )


                tuple(
                    sample_id,
                    fasta
                )
            }
            .ifEmpty {
                error """
                Samplesheet contains no samples:

                    ${params.input}
                """.stripIndent()
            }


            // Agrupa por sample_id para detectar identificadores repetidos antes de lanzar el pipeline
            .groupTuple()
            .map { sample_id, fastas ->

                if (fastas.size() != 1) {
                    error """
                    Duplicate sample_id in samplesheet:

                        ${sample_id}

                    Each sample_id must be unique.
                    """.stripIndent()
                }


                tuple(
                    sample_id,
                    fastas[0]
                )
            }
    }


    // Si recibe un único FASTA, lo convierte al mismo formato usado por el samplesheet
    else {

        if (!params.sample_id) {
            error """
            Missing required parameter: --sample_id
            When using:
                --fasta genome.fna
            you must also provide:
                --sample_id sample_name
            Example:
                nextflow run main.nf \\
                    --fasta genome.fna \\
                    --sample_id ecoli_test \\
                    -profile local
            """.stripIndent()
        }


        def sample_id = params.sample_id
            .toString()
            .trim()


        if (
            !sample_id.matches(
                '[A-Za-z0-9][A-Za-z0-9_.-]*'
            )
        ) {
            error """
            Invalid sample_id: ${sample_id}

            Allowed characters:

                A-Z a-z 0-9 _ - .

            The sample_id must start with a letter or number.
            """.stripIndent()
        }


        def fasta = file(
            params.fasta,
            checkIfExists: true
        )


        // Crea un canal con la misma estructura tuple(sample_id, fasta)
        samples_ch = channel.of(
            tuple(
                sample_id,
                fasta
            )
        )
    }


    // A partir de este punto ambos modos de entrada son iguales: samples_ch siempre emite tuple(sample_id, fasta)


    // Convierte la ruta de la base de datos de Bakta en un valor reutilizable por todas las muestras
    bakta_db_ch = channel.value(
        params.bakta_db
    )


    // Ejecuta Bakta sobre cada genoma y obtiene el JSON estructural y el FASTA de proteínas
    bakta = BAKTA(
        samples_ch,
        bakta_db_ch
    )


    // El JSON alimentará el parser de Bakta y el FAA será la entrada de las herramientas funcionales
    bakta_json_ch = bakta.results.map {
        sample_id,
        json,
        _faa ->

        tuple(
            sample_id,
            json
        )
    }


    // Se queda solo con el FASTA de proteínas de Bakta, que será la entrada de las herramientas funcionales
    bakta_faa_ch = bakta.results.map {
        sample_id,
        _json,
        faa ->

        tuple(
            sample_id,
            faa
        )
    }


    // Convierte el JSON de Bakta en las tablas canónicas de features, annotations y evidence
    parsed = PARSE_BAKTA(
        bakta_json_ch
    )


    // Comparte la base de datos de eggNOG entre todas las muestras sin copiarla dentro de cada tarea
    eggnog_db_ch = channel.value(
        params.eggnog_db
    )


    // Comparte la base de datos de KofamScan entre todas las muestras
    kofam_db_ch = channel.value(
        params.kofam_db
    )


    // Comparte el catálogo local utilizado para las relaciones KO-EC
    kofam_ko_list_ch = channel.value(
        file(
            "${params.kofam_db}/ko_list",
            checkIfExists: true
        )
    )


    // Comparte la ruta de la instalación externa de InterProScan entre todas las muestras
    interproscan_path_ch = channel.value(
        params.interproscan_path
    )


    // Carga la configuración de PathoFact que forma parte del propio repositorio
    pathofact_config_ch = channel.value(
        file(
            "${projectDir}/conf/pathofact.yaml",
            checkIfExists: true
        )
    )


    // Carga el script que adapta y prepara PathoFact para ejecutarlo dentro de BactAnnot
    pathofact_prepare_script_ch = channel.value(
        file(
            "${projectDir}/scripts/prepare_pathofact_runtime.py",
            checkIfExists: true
        )
    )

    // Carga los archivos auxiliares que necesita la capa de compatibilidad con PathoFact
    pathofact_assets_ch = channel.value(
        file(
            "${projectDir}/assets/pathofact",
            checkIfExists: true
        )
    )


    // Comparte la ruta de la instalación externa de PathoFact
    pathofact_path_ch = channel.value(
        params.pathofact_path
    )

    // Comparte la ruta de las bases de datos externas de PathoFact
    pathofact_db_ch = channel.value(
        params.pathofact_db
    )

    // Comparte la ubicación de los entornos Conda que utiliza PathoFact
    pathofact_conda_prefix_ch = channel.value(
        params.pathofact_conda_prefix
    )


    // Ejecuta eggNOG sobre las proteínas de Bakta para obtener anotaciones por ortología
    eggnog = EGGNOG(
        bakta_faa_ch,
        eggnog_db_ch
    )


    // Analiza las proteínas de Bakta buscando dominios, firmas funcionales y términos GO
    interproscan = INTERPROSCAN(
        bakta_faa_ch,
        interproscan_path_ch
    )


    // Busca asignaciones KEGG Orthology en las proteínas de Bakta
    kofamscan = KOFAMSCAN(
        bakta_faa_ch,
        kofam_db_ch
    )


    // Analiza las proteínas de Bakta para obtener información de virulencia, toxinas y AMR
    pathofact = PATHOFACT(
        bakta_faa_ch,
        pathofact_config_ch,
        pathofact_prepare_script_ch,
        pathofact_assets_ch,
        pathofact_path_ch,
        pathofact_db_ch,
        pathofact_conda_prefix_ch
    )


    // El parser de eggNOG necesita las anotaciones y los seed orthologs, pero no el archivo de hits
    eggnog_for_parse_ch = eggnog.results.map {
        sample_id,
        annotations,
        seed_orthologs,
        _hits ->

        tuple(
            sample_id,
            annotations,
            seed_orthologs
        )
    }


    // Une la salida de eggNOG con las features de Bakta usando sample_id para recuperar las feature_uid correctas
    eggnog_parse_ch = eggnog_for_parse_ch
        .join(
            parsed.features
        )


    // Normaliza los resultados de eggNOG al formato común de annotations y evidence
    parsed_eggnog = PARSE_EGGNOG(
        eggnog_parse_ch
    )


    // Une los resultados de InterProScan con las features de Bakta de la misma muestra
    interproscan_parse_ch = interproscan.results
        .join(
            parsed.features
        )


    // Convierte los resultados de InterProScan a las tablas canónicas de BactAnnot
    parsed_interproscan = PARSE_INTERPROSCAN(
        interproscan_parse_ch
    )


    // Une los resultados de KofamScan con las features de Bakta de la misma muestra
    kofamscan_parse_ch = kofamscan.results
        .join(
            parsed.features
        )


    // Convierte los resultados de KofamProScan a las tablas canónicas de BactAnnot
    parsed_kofamscan = PARSE_KOFAMSCAN(
        kofamscan_parse_ch
    )


    // Une todos los resultados de PathoFact con las features de Bakta de la misma muestra
    pathofact_parse_ch = pathofact.results
        .join(
            parsed.features
        )


    // Convierte los resultados de PathoFact a las tablas canónicas de BactAnnot
    parsed_pathofact = PARSE_PATHOFACT(
        pathofact_parse_ch
    )


    // Junta las tablas de anotaciones de todas las herramientas y las agrupa por muestra
    normalized_annotations_ch = parsed.annotations
        .mix(
            parsed_eggnog.annotations,
            parsed_interproscan.annotations,
            parsed_pathofact.annotations,
            parsed_kofamscan.annotations
        )
        .groupTuple()


    // Añade la tabla de features para que la integración pueda validar cada anotación
    annotation_integration_ch = normalized_annotations_ch
        .join(
            parsed.features
        )


    // Integra todas las anotaciones sin perder la herramienta de la que procede cada una
    annotation_integration = INTEGRATE_ANNOTATIONS(
        annotation_integration_ch
    )


    // Junta las tablas de evidencia de todas las herramientas y las agrupa por muestra
    normalized_evidence_ch = parsed.evidence
        .mix(
            parsed_eggnog.evidence,
            parsed_interproscan.evidence,
            parsed_pathofact.evidence,
            parsed_kofamscan.evidence
        )
        .groupTuple()


    // Añade las features de Bakta para validar a qué elemento pertenece cada evidencia
    evidence_integration_ch = normalized_evidence_ch
        .join(
            parsed.features
        )


    // Integra toda la evidencia conservando su fuente original
    evidence_integration = INTEGRATE_EVIDENCE(
        evidence_integration_ch
    )


    // Reúne anotaciones, evidencias y features de cada muestra para construir el consenso
    consensus_input_ch = annotation_integration.annotations
        .join(
            evidence_integration.evidence
        )
        .join(
            parsed.features
        )


    // Construye claims por namespace y relaciones funcionales documentadas
    consensus_build = BUILD_CONSENSUS(
        consensus_input_ch,
        kofam_ko_list_ch
    )


    // Concordance solo necesita las anotaciones integradas y las features de Bakta
    concordance_input_ch = annotation_integration.annotations
        .join(
            parsed.features
        )


    // Compara por pares qué identificadores funcionales aporta cada herramienta
    concordance_build = BUILD_CONCORDANCE(
        concordance_input_ch
    )

    // El informe recibe archivos de los canales, antes de su publicación en results.
    report_output_ch = channel.empty()

    if (!params.skip_report) {

        concordance_report_ch = concordance_build.concordance.map {
            sample_id, concordance_table, _summary ->
            tuple(sample_id, concordance_table)
        }

        report_input_ch = parsed.features
            .join(annotation_integration.annotations, failOnDuplicate: true, failOnMismatch: true)
            .join(evidence_integration.evidence, failOnDuplicate: true, failOnMismatch: true)
            .join(consensus_build.consensus, failOnDuplicate: true, failOnMismatch: true)
            .join(concordance_report_ch, failOnDuplicate: true, failOnMismatch: true)

        report_build = BUILD_REPORT(
            report_input_ch,
            channel.value(file("${projectDir}/report/bactannot_report.qmd", checkIfExists: true)),
            channel.value(file("${projectDir}/src", checkIfExists: true))
        )

        report_output_ch = report_build.reports
    }


    publish:

    // Expone los resultados que queremos conservar al terminar el workflow

    // Bakta
    bakta_raw = bakta.results
    features = parsed.features
    annotations = parsed.annotations
    evidence = parsed.evidence


    // eggNOG
    eggnog_raw = eggnog.results
    eggnog_annotations = parsed_eggnog.annotations
    eggnog_evidence = parsed_eggnog.evidence


    // InterProScan
    interproscan_raw = interproscan.results
    interproscan_annotations = parsed_interproscan.annotations
    interproscan_evidence = parsed_interproscan.evidence


    // KofamScan
    kofamscan_raw = kofamscan.results
    kofamscan_annotations = parsed_kofamscan.annotations
    kofamscan_evidence = parsed_kofamscan.evidence


    // PathoFact
    pathofact_raw = pathofact.results
    pathofact_annotations = parsed_pathofact.annotations
    pathofact_evidence = parsed_pathofact.evidence


    // Integración
    integrated_annotations = annotation_integration.annotations
    integrated_evidence = evidence_integration.evidence


    // Consensus
    consensus = consensus_build.consensus


    // Concordance
    concordance = concordance_build.concordance

    reports = report_output_ch
}


output {

    // Define en qué carpeta de results se guarda cada tipo de salida

    // Bakta
    bakta_raw {
        path 'bakta/raw'
    }

    features {
        path 'bakta/features'
    }

    annotations {
        path 'bakta/annotations'
    }

    evidence {
        path 'bakta/evidence'
    }


    // eggNOG
    eggnog_raw {
        path 'eggnog/raw'
    }

    eggnog_annotations {
        path 'eggnog/annotations'
    }

    eggnog_evidence {
        path 'eggnog/evidence'
    }


    // InterProScan
    interproscan_raw {
        path 'interproscan/raw'
    }

    interproscan_annotations {
        path 'interproscan/annotations'
    }

    interproscan_evidence {
        path 'interproscan/evidence'
    }


    // KofamScan
    kofamscan_raw {
        path 'kofamscan/raw'
    }

    kofamscan_annotations {
        path 'kofamscan/annotations'
    }

    kofamscan_evidence {
        path 'kofamscan/evidence'
    }


    // PathoFact
    pathofact_raw {
        path 'pathofact/raw'
    }

    pathofact_annotations {
        path 'pathofact/annotations'
    }

    pathofact_evidence {
        path 'pathofact/evidence'
    }


    // Integración
    integrated_annotations {
        path 'integration/annotations'
    }

    integrated_evidence {
        path 'integration/evidence'
    }


    // Consensus
    consensus {
        path 'consensus'
    }


    // Concordance
    concordance {
        path 'concordance'
    }

    reports {
        path 'reports'
        mode 'copy'
        overwrite true
    }
}