import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

val in = new BufferedReader(new FileReader(new File("/home/aeonsim/refs/DamonaGenos20160407.txt")))
//val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(new File("/home/projects/bos_taurus/damona/vcfs/GATK-HC-3.4-Damona_12_Jan_2016_sorted.vcf.gz")))))
val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(new File("/home/aeonsim/refs/BBB-seqdHDs.AC1plus.vcf.gz")))))

val out = new BufferedWriter(new FileWriter(new File("/home/aeonsim/Valdiation_50K.vcf")))

val bases = new HashMap[String,Tuple2[String,String]]

var line = in_vcf.readLine.split("\t")

while(in_vcf.ready){
	
	while (line(0)(0) == '#') line = in_vcf.readLine.split("\t")
	println("adding " + s"${line(0)}:${line(1)}")
	bases += s"${line(0)}:${line(1)}" -> Tuple2(line(3),line(4))
	line = in_vcf.readLine.split("\t")
}
in_vcf.close


out.write("##fileformat=VCFv4.1\n")
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

val hetStates = List("RA","AR")
var cline = in.readLine.split("\t")

for (i <- 4 to cline.size -1){
	out.write("\t" + cline(i))
}

out.write("\n")

while (in.ready){
	cline = in.readLine.split("\t")
	val pos = s"chr${cline(2)}:${cline(3)}}"
	if (bases.contains(pos)){
		val cbases = bases(pos)
		out.write(s"chr${cline(2)}\t${cline(3)}\t${cline(0)}\t${cbases._1}\t${cbases._2}\t.\t.\tGT")

		for (i <- 4 to cline.size -1){
			if (hetStates.contains(cline(i))) out.write("\t0/1") else if (cline(i) == "00") out.write("\t./.") else if (cline(i) == "AA" ) out.write("\t1/1") else if (cline(i) == "RR") out.write("\t0/0") else out.write("\t.")
			out.write("\n")
		}
	}

}

out.close