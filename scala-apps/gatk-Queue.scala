import java.io.File
import scala.util.Random

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction

class BAM2VCF extends QScript {
  qscript =>

  /****************************************************************************
   * Required Parameters
   *****************************************************************************/

  @Input(doc="The reference file for the bam files.", shortName="R", fullName="reference_sequence")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to relalign", shortName="I", fullName="input_file")
  var myBam: File = _

  /****************************************************************************
   * Optional Parameters
   *****************************************************************************/

  @Argument(doc="scatter parameter", shortName="P", fullName="scatter_parameter", required=false)
  var scatter:  Int = _

  @Argument(doc="nt parameter", shortName="N", fullName="num_threads", required=false)
  var nt:  Int = _

  @Argument(doc="nct parameter", shortName="C", fullName="num_cpu_threads_per_data_thread", required=false)
  var nct:  Int = _

  @Input(doc="Intervals to realign", shortName="L", required = false)
  var intervals: File = _

  @Input(doc="dbSNP file", shortName="D", fullName="dbsnp", required=false)
  var dbsnp: File = _

  @Input(doc="Known in/del VCF file", shortName="known", required=false)
  var known: File = _

  @Input(doc="BQSR recalibration table", shortName="BQSR", required=false)
  var BQSR: File = _

  @Input(doc="First recalibration table", shortName="before", required=false)
  var before: File = _

  @Input(doc="Second recalibration table; BQSR with first recalibration table", shortName="after", required=false)
  var after: File = _

  /****************************************************************************
   * CommonArguments
   *****************************************************************************/

  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 20
  }


  /****************************************************************************
   * Main script
   *****************************************************************************/

  def script() {

    val targetCreator = new RealignerTargetCreator with CommonArguments
    val indelRealigner = new IndelRealigner with CommonArguments
    val bqsr = new BaseRecalibrator with CommonArguments
    val bqsr2 = new BaseRecalibrator with CommonArguments
    val analyzeCovariates = new AnalyzeCovariates with CommonArguments
    val applyRecalibration = new PrintReads with CommonArguments
    val hc = new HaplotypeCaller with CommonArguments
	val combine = new CombineGVCFs with CommonArguments
	val genotype = new GenotypeGVCFs with CommonArguments
	
	
    targetCreator.input_file +:= qscript.myBam
    targetCreator.out = swapExt(myBam, ".bam", ".intervals")
    targetCreator.nt = qscript.nt
    targetCreator.known +:= qscript.known

    indelRealigner.input_file +:= qscript.myBam
    indelRealigner.targetIntervals = targetCreator.out
    indelRealigner.out = swapExt(myBam, ".dd.bam", ".dd.ir.bam")
    indelRealigner.known +:= qscript.known

    bqsr.input_file +:= indelRealigner.out
    bqsr.out = swapExt(indelRealigner.out, ".dd.ir.bam", ".bqrecal")
    bqsr.knownSites +:= new File("/scratch/ulg/genan/rgularte/mus_musculus/Ensembl/Grcm38/Annotation/Variation/Mus_musculus.sorted.vcf") //Seq(qscript.KS)
    bqsr.nct = qscript.nct

    bqsr2.input_file +:= indelRealigner.out
    bqsr2.out = swapExt(bqsr.out, ".bqrecal", ".bqrecal2")
    bqsr2.knownSites +:= new File("/scratch/ulg/genan/rgularte/mus_musculus/Ensembl/Grcm38/Annotation/Variation/Mus_musculus.sorted.vcf") //Seq(qscript.KS)
    bqsr2.BQSR = bqsr.out
    bqsr2.nct = qscript.nct

    analyzeCovariates.before = bqsr.out
    analyzeCovariates.after = bqsr2.out
    analyzeCovariates.csv = swapExt(bqsr.out, ".dd.ir.bam", ".recal.csv")
    analyzeCovariates.plots = swapExt(bqsr.out, ".dd.ir.bam", ".bqsr.pdf")

    applyRecalibration.input_file +:= indelRealigner.out
    applyRecalibration.BQSR = bqsr2.out
    applyRecalibration.out = swapExt(indelRealigner.out, ".dd.ir.bam", ".dd.ir.bqsr.bam")
    applyRecalibration.nct = qscript.nct

    hc.scatterCount = qscript.scatter
    hc.input_file +:= applyRecalibration.out
    hc.nct = qscript.nct
    hc.pairHMM = org.broadinstitute.gatk.utils.pairhmm.PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING
    hc.pcr_indel_model = org.broadinstitute.gatk.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE
//    hc.out = swapExt(qscript.myBam, ".dd.ir.bqsr.bam", ".vcf")
    hc.out = swapExt(qscript.myBam, ".dd.ir.bqsr.bam", ".g.vcf")
    hc.emitRefConfidence = org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
    hc.variant_index_type = org.broadinstitute.gatk.utils.variant.GATKVCFIndexType.LINEAR
    hc.variant_index_parameter = 128000

    add(targetCreator, indelRealigner, bqsr, bqsr2, analyzeCovariates, applyRecalibration, hc)

  }

}

