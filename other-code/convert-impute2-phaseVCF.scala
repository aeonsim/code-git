import java.io._

val anmls = scala.io.Source.fromFile("test.effective-pop-300.out.sample").getLines.toArray
val data = new BufferedReader(new FileReader("test.effective-pop-300.out.haps"))

val out = new BufferedWriter(new FileWriter("impute2.vcf"))

out.write("##fileformat=VCFv4.1\n")

val names = new Array[String](anmls.size - 2 )
for (x <- 2 to (anmls.size -1)){
names(x-2) = anmls(x).split(" ")(1)
}

out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + names.reduceLeft[String]((a,b) => a + "\t" + b) + "\n")
while (data.ready){
val line = data.readLine.split(" ")
if (line(3) != line(4)){
out.write("Chr1\t" + line(2) + "\t.\t" + line(3) + "\t"+ line(4) + "\t.\t.\t.\tGT")
var pos = 5
while (pos <= (line.size -1)){
out.write( "\t" + line(pos) + "|" + line(pos + 1))
pos += 2
}
out.write("\n")
}
}
out.close