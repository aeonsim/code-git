import java.io._

val input = new BufferedReader(new FileReader("AGPAT6-pop-UG-preprocessed.vcf"))

var line = input.readLine
while (line(1) == '#'){
line = input.readLine
}

var cline = line.split("\t")
val numAn = cline.size - 9
//0 Chrom 8 FORMAT 9 -> Animal ID
var animals = new Array[List[Int]](numAn)
var names = new Array[String](numAn)

for (cAn <- 0 to (numAn -1)){
animals(cAn) = Nil
names(cAn) = cline(cAn + 9)
}

//ID GT location
cline = input.readLine.split("\t")
val dpPos = cline(8).split(":").indexOf("DP")

while (input.ready){
for (cAn <- (0 to (numAn - 1))){
val geno = cline(cAn + 9).split(":")
if (geno.size >= 2){
animals(cAn) = geno.apply(dpPos).toInt :: animals(cAn)
}
}
cline = input.readLine.split("\t")
}


for (cAn <- (0 to (numAn -1))){
println(names(cAn) + ": " + animals(cAn).sum / (if(animals(cAn).size > 0){animals(cAn).size}else{1}))
}