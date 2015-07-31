object checker {

def main(args: Array[String]): Unit = {

import java.io._
import net.sf.samtools.util._

if (args.size != 2){
println("checker VCFfile.vcg.gz PEDfile.ped")
System.exit(1)
}

val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(args(0)))))
val in_ped = new BufferedReader(new FileReader(args(1)))

var ped : List[String] = Nil

while (in_ped.ready){
val line = in_ped.readLine.split("\t")
ped = line(1) :: ped
}

var vcfLine = in_vcf.readLine.split("\t")

while (vcfLine(0)(0) == '#' && vcfLine(0)(1) != 'C') vcfLine = in_vcf.readLine.split("\t")

for (i <- vcfLine){
if (! ped.contains(i)){println(i + " Not in Pedigree File")}
}


}//Main
}
