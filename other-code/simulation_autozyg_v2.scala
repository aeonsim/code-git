import java.io._
import scala.collection.mutable.HashMap
import org.apache.commons.math3.distribution._

object gens {

val rand = scala.util.Random

def main(args:Array[String]) : Unit = {

if (args.size != 5){
println("gens pedIn vcfIn pedOut vcfOut matings")
scala.sys.exit
}

val gensIn = new BufferedReader(new FileReader(args(0)))
val vcfIn = new BufferedReader(new FileReader(args(1)))

val gensOut = new BufferedWriter(new FileWriter(args(2)))
val vcfOut = new BufferedWriter(new FileWriter(args(3)))
val matingsOut = new BufferedWriter(new FileWriter(args(4)))


/*/////////////////////
val gensIn = new BufferedReader(new FileReader("ped.txt"))
val gensOut = new BufferedWriter(new FileWriter("gens1.txt"))
val vcfIn = new BufferedReader(new FileReader("g0.vcf"))
val vcfOut = new BufferedWriter(new FileWriter("vcfout.vcf"))
*/

// ID, Sex, Generations
var ped: List[Tuple3[Int,String,Int]] = Nil

while (gensIn.ready){
val cline = gensIn.readLine.split("\t")
ped = (cline(0).toInt,cline(1),(cline(2).toInt + 1)) :: ped
}

val males = ped.filter(a => a._2 == "M" && a._3 <= 3)
val females = ped.filter(a => a._2 == "F" && a._3 <= 3)

/*
* Works out number of Children Per Sire
*/

val sireChild = females.size / males.size
val maxAnID = ped.map(_._1).max + females.size
var curAnml = ped.map(_._1).max + 1
var sons = 0

while (curAnml <= maxAnID){
if (sons < 2) {
ped = (curAnml,"M",0) :: ped
sons += 1
} else {
ped = (curAnml,"F",0) :: ped
}
curAnml += 1
}

for(an <- ped.reverse){
gensOut.write(an._1 + "\t" + an._2 + "\t" + an._3 + "\n")
}
gensOut.close


/*
* Mut and CO Code
*/

val chrsIn = new BufferedReader(new FileReader("chrs.txt"))

val mutRate = 1.2E-8
val co = 0.07 //per 10Mb
val cobase = 0.48
val Chroms = new Array[Tuple2[String,Int]] (29)
val probs = new Array[PoissonDistribution] (29)
val coprobs = new Array[PoissonDistribution] (29)

for (num <- 0 to 28) {
val cline = chrsIn.readLine.split("\t")
Chroms(num) = (cline(0),cline(1).toInt) 
probs(num) = new PoissonDistribution(mutRate*Chroms(num)._2)
coprobs(num) = new PoissonDistribution((Chroms(num)._2/10000000)* co + cobase)
}
chrsIn.close

/* 
* Denovo Mutation Generator, Creates Denovo mutations spread across all Chrs
* Base on Chr size and average Mut rate, # of Mutations per individual
* follow poisson dist with the mean number of mutations for that size chromosome
*/
val rand = new scala.util.Random

//Chr, Pos, Anml
var mutations : List[Tuple3[Int,Int,Int]] = Nil
for(dams <- females){
var chr = 0 
while (chr < 29 ){
val numMut = probs(chr).sample
var muts = 0
while (muts < numMut){
//Chromosome Name, Position, Animal ID
mutations = (chr,rand.nextInt(Chroms(chr)._2 + 1),dams._1) :: mutations
muts += 1
}
chr += 1
}
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
for(dams <- females){
var chr = 0 
while (chr < 29 ){
val numCO = coprobs(chr).sample
var comuts = 0
while (comuts < numCO){
//Chromosome Name
coEvents = (chr,rand.nextInt(Chroms(chr)._2 + 1),dams._1) :: coEvents
comuts += 1
}
chr += 1
}
}

//var muts = Array[Tuple3[String, Int, Int]]
/*
*/



/*
* Mattings List
*/
//val matings = scala.util.Random.shuffle(females).grouped(sireChild).toArray

var vcfLine = vcfIn.readLine.split("\t")

//Read through headers and write them out
while ((vcfLine(0)(0) == '#')&&(vcfLine(0)(1) == '#')){
vcfOut.write(vcfLine.reduceLeft[String]((a,b) => a + "\t" + b) + "\n")
vcfLine = vcfIn.readLine.split("\t")
}

//Write out #Chrom Line with new animal ID's
vcfOut.write(vcfLine(0) + "\t" + vcfLine(1)+ "\t" +vcfLine(2)+ "\t" +vcfLine(3)+ "\t" +vcfLine(4)+ "\t" +vcfLine(5)+ "\t" +vcfLine(6)+ "\t" +vcfLine(7) + "\t" + vcfLine(8) + "\t")
vcfOut.write(ped.reverse.map(_._1.toString).reduceLeft[String]((a,b) => a + "\t" + b) + "\n")

/*
* Select Phase to Use for Each animal
*/

var femalePhase = new HashMap[Int,Int]
var malePhase = new HashMap[Int,Int]
for (cur <- females){
femalePhase += cur._1 -> rand.nextInt(2)
}
for (cur <- males){
malePhase += cur._1 -> rand.nextInt(2)
}
var sire4dam = new HashMap[Int,Int]
var curChild = 0
for (cur <- females){
curChild += 1
sire4dam += cur._1 -> males(rand.nextInt(males.size))._1
matingsOut.write((vcfLine.size + curChild) + "\t" + sire4dam(cur._1) + "\t" + cur._1 + "\n")
}
matingsOut.close

/*
* Setup Mutation and CO points
*/
vcfLine = vcfIn.readLine.split("\t")
var curChr = vcfLine(0)
var chrIndx = vcfLine(0).toUpperCase.split("CHR")(1).toInt - 1
var cmuts = mutations.filter(_._1 == chrIndx).sortBy(_._2).toArray
var workingMut = 0
var comuts = coEvents.filter(_._1 == chrIndx).sortBy(_._2).toArray
var workingCO = 0

while(vcfIn.ready && (vcfLine(0).toUpperCase.split("CHR")(1) != "X")){

if(curChr != vcfLine(0)){
curChr = vcfLine(0)
chrIndx = vcfLine(0).toUpperCase.split("CHR")(1).toInt - 1
cmuts = mutations.filter(_._1 == chrIndx).sortBy(_._2).toArray
comuts = coEvents.filter(_._1 == chrIndx).sortBy(_._2).toArray
workingMut = 0
workingCO = 0
}

if((workingMut < cmuts.size) && (vcfLine(1).toInt > cmuts(workingMut)._2)){
println("Denovo " +cmuts(workingMut)._2)
vcfOut.write(vcfLine(0) +"\t"+ cmuts(workingMut)._2 + "\tDENOVO\tA\tT\t.\tPASS\t.\tGT")
for (num <- 0 to ped.size){
if (cmuts(workingMut)._3 == num) vcfOut.write("\t0/1") else vcfOut.write("\t0/0")
}
vcfOut.write("\n")
workingMut += 1
}

vcfOut.write(vcfLine.reduceLeft[String]((a,b) => a + "\t" + b))

for (dam <- females.reverse){
curChild += 1
//CO controls
val damPro = vcfLine(dam._1 + 8)
val sirePro = vcfLine(sire4dam(dam._1) +8)

if((workingCO < comuts.size) && (vcfLine(1).toInt > comuts(workingCO)._2)){
	println("CO Event " + comuts(workingCO)._2)
	if (dam._1 == comuts(workingCO)._3){
		if (femalePhase(dam._1) == 0) femalePhase(dam._1) = 1
		else femalePhase(dam._1) = 0
	}
	if (sire4dam(dam._1) == comuts(workingCO)._3){
		if (malePhase(sire4dam(dam._1)) == 0) malePhase(dam._1) = 1
		else malePhase(dam._1) = 0
	}
	workingCO += 1
}
val childPro = damPro.split("/")(femalePhase(dam._1)) + "/" + sirePro.split("/")(malePhase(sire4dam(dam._1)))
vcfOut.write("\t" + childPro)
}
vcfOut.write("\n")
vcfLine = vcfIn.readLine.split("\t")
}
vcfIn.close
vcfOut.close
}
//end def main
}

//end object