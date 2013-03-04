import scala.io.Source._
import java.io._

object prepVcf {

def main ( args: Array[String]){

val vcf = new BufferedReader(new FileReader(args(0)))

//Load VCF
//val in = fromFile("genos.txt").getLines.toArray
val out = new BufferedWriter(new FileWriter("homozygosity_id.txt"))
val mafs = new BufferedWriter(new FileWriter("mafs.txt"))
val indx = new BufferedWriter(new FileWriter("snp_indx.txt"))
val animals = new BufferedWriter(new FileWriter("anml_indx.txt"))

var cnt = 0

var cl = vcf.readLine.split("\t")
while (cl(0).apply(0) == '#') {
if (cl(0).apply(1) == 'C'){
for (anmls <- 9 to (cl.size -1)){
animals.write((anmls - 8) + "\t" + cl(anmls) + "\n")
}
animals.close
}
cl = vcf.readLine.split("\t")
}

val GTpos = cl(8).split(":").indexOf("GT")

while (vcf.ready){
//Write index Line
indx.write((cnt + 1) + "\t" + cl(0) + "\t" + cl(1) + "\n")
//Write Homzygosity map
out.write((cnt+1) + "")
var AA, BB = 0
for (anmls <- 9 to (cl.size -1)){
cl(anmls).split(":").apply(GTpos) match {
 case "0/0" => out.write(" 1.0"); AA += 2
 case "0/1" => out.write(" 0.0"); AA += 1; BB += 1
 case "1/0" => out.write(" 0.0"); AA += 1; BB += 1
 case "1/1" => out.write(" 1.0"); BB += 2
 case "./." => out.write(" 2.0")
 case _ => out.write(" 2.0")
 }
} //for end
out.write("\n")
 val MAF = if ((AA + BB) == 0) 0 else AA.toFloat / (AA + BB).toFloat
mafs.write((cnt + 1) + "  " + MAF + "\n")
 cnt += 1
 cl = vcf.readLine.split("\t")
}//while end
vcf.close
out.close
mafs.close
indx.close
}
}