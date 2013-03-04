import org.apache.commons.math3.distribution._
import java.io._
import scala.collection.mutable.HashMap

object mutnCO {

def main(args: Array[String]) : Unit = {

if (args.size != 2){
println("mutnCO in.vcf out.vcf")
scala.sys.exit
}
//File format Chromosome, Size
val chrsIn = new BufferedReader(new FileReader("chrs.txt"))

val mutRate = 1.2E-8
val co = 0.07 //per 10Mb
val cobase = 0.48
val Chroms = new Array[Tuple2[String,Int]] (29)
val probs = new Array[PoissonDistribution] (29)
val coprobs = new Array[PoissonDistribution] (29)

val out = new BufferedWriter(new FileWriter(args(1)))

for (num <- 0 to 28) {
val cline = chrsIn.readLine.split("\t")
Chroms(num) = (cline(0),cline(1).toInt) 
probs(num) = new PoissonDistribution(mutRate*Chroms(num)._2)
coprobs(num) = new PoissonDistribution((Chroms(num)._2/10000000)* co + cobase)
}
chrsIn.close

/* 
* Loads Gender data for determining matings and snps from VCF for generating profile
*/

//val gender =  new BufferedReader(new FileReader("gender.txt"))
val vcf = new BufferedReader(new FileReader(args(0)))
//val snpsIn = new BufferedReader(new FileReader("snp_indx.txt"))

//val Index = new HashMap[String, Int]
/*
while (snpsIn.ready){
val cline = snpsIn.readLine.split("\t")
Index += (cline(1) + ":" + cline(2)) -> cline(0).toInt
}
*/
//snpsIn.close

//Deal with VCFs
/*
* Load VCF file and determine Position of GT's in the Animal details
*/

val animals = new BufferedWriter(new FileWriter("anml_indx.txt"))
var cl = vcf.readLine.split("\t")
while (cl(0).apply(0) == '#') {
out.write(cl.reduceLeft[String]((a,b)=> a +  "\t" + b) + "\n")
if (cl(0).apply(1) == 'C'){
for (anmls <- 9 to (cl.size -1)){
animals.write((anmls - 8) + "\t" + cl(anmls) + "\n")
}
animals.close
}
cl = vcf.readLine.split("\t")
}
val GTpos = cl(8).split(":").indexOf("GT")

/*
* Matings generator, selects males and randomly mates them with multiple
* females. Mattings output in a text file.
*
*/

/* 
* Denovo Mutation Generator, Creates Denovo mutations spread across all Chrs
* Base on Chr size and average Mut rate, # of Mutations per individual
* follow poisson dist with the mean number of mutations for that size chromosome
*/
val numAnmls = cl.size - 9
val rand = new scala.util.Random

//Chr, Pos, Anml
var mutations : List[Tuple3[Int,Int,Int]] = Nil
var anmls = 0
while(anmls < numAnmls){
var chr = 0 
while (chr < 29 ){
val numMut = probs(chr).sample
var muts = 0
while (muts < numMut){
//Chromosome Name
mutations = (chr,rand.nextInt(Chroms(chr)._2 + 1),anmls) :: mutations
muts += 1
}
chr += 1
}
anmls += 1
}

/*
* Recombination event Gen, similar to Denovo gen
* Calc average number of crossovers per chromosome and implement.
* Treating 1/1 as phased for simplification. Using Poisson Dist based
* on per Chromosome averages from:
*"Genetic variants in REC8, RNF212 and PRDM9 influence male recombination in cattle."
* Cynthia Sandor#1, Wanbo Li1, Wouter Coppieters, Tom Druet, Carole Charlier & Michel Georges.
*/

var coEvents : List[Tuple3[Int,Int,Int]] = Nil
var coanmls = 0
while(coanmls < numAnmls){
var chr = 0 
while (chr < 29 ){
val numCO = coprobs(chr).sample
var comuts = 0
while (comuts < numCO){
//Chromosome Name
coEvents = (chr,rand.nextInt(Chroms(chr)._2 + 1),coanmls) :: coEvents
comuts += 1
}
chr += 1
}
coanmls += 1
}

//var muts = Array[Tuple3[String, Int, Int]]
var curChr = cl(0)
var chrIndx = cl(0).toUpperCase.split("CHR")(1).toInt - 1
var cmuts = mutations.filter(_._1 == chrIndx).sortBy(_._2).toArray
var workingMut = 0
var comuts = coEvents.filter(_._1 == chrIndx).sortBy(_._2).toArray
var workingCO = 0

/*
* Outputs the updated VCF containing Denovo Mutations and cross over events.
*/

/*
* Outputs the updated VCF containing Denovo Mutations and cross over events.
*/

var anmlCOstate = new Array[Boolean](numAnmls)

println(numAnmls)

def writeLine(): Unit = {
out.write(cl(0) + "\t" + cl(1)+ "\t" +cl(2)+ "\t" +cl(3)+ "\t" +cl(4)+ "\t" +cl(5)+ "\t" +cl(6)+ "\t" +cl(7) + "\tGT")
for (coNUM <- 9 to (8 + numAnmls)){
val gts = cl(coNUM).split(":").apply(GTpos).split("/")
if (anmlCOstate(coNUM - 9)){
out.write("\t" + gts(1) + "/" + gts(0))
} else {
out.write("\t" + gts(0) + "/" + gts(1))
}
}
out.write("\n")
}


while (vcf.ready && (cl(0).toUpperCase.split("CHR")(1) != "X")){
if(curChr != cl(0)){
curChr = cl(0)
chrIndx = cl(0).toUpperCase.split("CHR")(1).toInt - 1
cmuts = mutations.filter(_._1 == chrIndx).sortBy(_._2).toArray
comuts = coEvents.filter(_._1 == chrIndx).sortBy(_._2).toArray
anmlCOstate = new Array[Boolean](numAnmls)
workingMut = 0
workingCO = 0
}

if (workingCO < comuts.size){
if (cl(1).toInt > comuts(workingCO)._2) {
anmlCOstate(comuts(workingCO)._3) = anmlCOstate(comuts(workingCO)._3).unary_!
workingCO += 1
}
}
if (workingMut == cmuts.size) {
writeLine
} else {
if (cl(1).toInt < cmuts(workingMut)._2){
writeLine
} else {
out.write(cl(0) +"\t"+ cmuts(workingMut)._2 + "\t.\tN\tN\t.\tPASS\t.\tGT")
for (num <- 0 to numAnmls){
if (cmuts(workingMut)._3 == num) out.write("\t0/1") else out.write("\t0/0")
}
out.write("\n")
writeLine
workingMut += 1
}
}
cl = vcf.readLine.split("\t")
}

out.close

}
}