import java.io._

val input = new BufferedReader(new FileReader("GATK-LIC-Het4-Hom4-DP1k-3k-FS-20.vcf"))
val out = new BufferedWriter(new FileWriter("Genos4m"))

var line = input.readLine
while (line(1) == '#'){
line = input.readLine
}

var cline = line.split("\t")
//0 Chrom 8 FORMAT 9 -> Animal ID
var animals = new Array[List[Int]](cline.size - 9)
var names = new Array[String](cline.size - 9)

for (cAn <- 9 to (cline.size - 1)){
animals(cAn - 9) = Nil
names(cAn - 9) = cline(cAn)
}

//ID GT location

cline = input.readLine.split("\t")
val gtPos = cline(8).split(":").indexOf("GT")

while (input.ready){
for (cAn <- 9 to (cline.size - 1)){
val gt = cline(cAn).split(":")(gtPos)
if ((gt.size == 3) && (gt(2) != '2')){
animals(cAn - 9) =  (gt(0).toString.toInt + 1) :: (gt(2).toString.toInt + 1) :: animals(cAn - 9)
} else {
animals(cAn - 9) = 0 :: 0 :: animals(cAn - 9)
}
}
cline = input.readLine.split("\t")
}

for (an <- 0 to (animals.size -1 )){
//out.write(names(an) + " ")
out.write(an + 1 + " ")
out.write(animals(an).reverse.map(_.toString).reduceLeft[String]((a,b) => a + " " + b))
out.write("\n")
}
println
out.close
