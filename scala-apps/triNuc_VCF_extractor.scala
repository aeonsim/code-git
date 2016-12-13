import htsjdk.samtools.reference.ReferenceSequence
import htsjdk.samtools.reference.ReferenceSequenceFile
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util._
import java.io._
import scala.collection.mutable.HashMap

val ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(args(0)))
val vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(1)))))

var vcfLine = vcf.readLine.split("\t")

var triNucCounts = new HashMap[String,Int]

while (vcfLine(0)(0) == '#') vcfLine = vcf.readLine.split("\t")

while (vcf.ready){
	val refBases = ref.getSubsequenceAt(vcfLine(0),vcfLine(1).toLong-1,vcfLine(1).toLong+1).getBases.map(_.toChar).mkString("").toUpperCase
	//println(refBases + " " + vcfLine(3) + " > " + vcfLine(4) )
	val altBases = if (vcfLine(3).size == 1 && vcfLine(4).size == 1 && vcfLine(3) == refBases(1).toString) refBases(0) + vcfLine(4) + refBases(2) else "INDEL"
	if(triNucCounts.contains(s"${refBases}>${altBases}")) triNucCounts(s"${refBases}>${altBases}") += 1 else triNucCounts += s"${refBases}>${altBases}" -> 1
	println(vcfLine(0) + ":" + vcfLine(1) + "\t" + refBases + ">" + altBases)
	vcfLine = vcf.readLine.split("\t")
}

ref.close
vcf.close

triNucCounts.foreach(s => println(s._1 + "\t" + s._2))