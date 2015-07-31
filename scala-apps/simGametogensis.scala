//JAVA_OPTS=-Xmx4g  scala -cp ~/Dropbox/PHD-WORK/pedigreeFilter/lib/commons-math3-3.3.jar

import java.io._

for (n <-  (1 to 100)){

val rand = scala.util.Random
val genomeSize = 2147000000
//val genomeSize = 2600000000L
val addSharedSoma = 8

/* Different differentiation events in embryo development */

val innerCellForm = 7
val innerCellPercent = 0.21
val ebiBlastFormation = 4
val epiBlastPercent = 0.5
val embEbpi = 4
val embEbpiPercent = 0.5


val pgcRepMale = 16
val pgcRepFemale = 18
val sprmgon = 500000
val oogon = 500000
val poi = new org.apache.commons.math3.distribution.PoissonDistribution(1)
val spermpoi = new org.apache.commons.math3.distribution.PoissonDistribution(2/23.0)
val mutFreqBlood = new scala.collection.mutable.HashMap[Int, Int]
val mutFreqFemaleGamete = new scala.collection.mutable.HashMap[Int, Int]
val mutFreqMaleGamete = new scala.collection.mutable.HashMap[Int, Int]
val yold = 5

def mutGen (num: Int): List[Int] = {
var count = 0
var result : List[Int] = Nil
while (count < num) {
result = rand.nextInt(genomeSize) :: result
count += 1
}
result
}

var cells: List[List[Int]] = List(mutGen(poi.sample))

for (i <- 1 to innerCellForm) cells = cells.map(e => mutGen(poi.sample) ::: e) ::: cells.map(e => mutGen(poi.sample) ::: e)

var newCells: List[List[Int]] = Nil

for (i <- 1 to (cells.size * innerCellPercent).toInt) newCells = cells(rand.nextInt(cells.size)) :: newCells

cells = newCells
newCells = Nil

for (i <- 1 to ebiBlastFormation) cells = cells.map(e => mutGen(poi.sample) ::: e) ::: cells.map(e => mutGen(poi.sample) ::: e)

for (i <- 1 to (cells.size * epiBlastPercent).toInt) newCells = cells(rand.nextInt(cells.size)) :: newCells

cells = newCells
newCells = Nil

for (i <- 1 to embEbpi) cells = cells.map(e => mutGen(poi.sample) ::: e) ::: cells.map(e => mutGen(poi.sample) ::: e)

for (i <- 1 to (cells.size * embEbpiPercent).toInt) newCells = cells(rand.nextInt(cells.size)) :: newCells

cells = newCells
newCells = Nil

var PGC: List[List[Int]] = Nil
for (i <- 1 to 6) PGC = cells(rand.nextInt(cells.size)) :: PGC
println(s"Selected ${PGC.size} PGCs, Soma is ${cells.size} cells")

for (i <- 1 to addSharedSoma) cells = cells.map(e => mutGen(poi.sample) ::: e) ::: cells.map(e => mutGen(poi.sample) ::: e)

for (i <- 1 to 34) PGC = cells(rand.nextInt(cells.size)) :: PGC
val blood = cells
cells = Nil
println(s"After ${addSharedSoma} additional reps have selected 34 additional PGCs, now have ${PGC.size} PGCs")

var PGCmA, PGCfA = Array(List(0))

for (i <- 1 to scala.math.max(pgcRepMale,pgcRepFemale)) {
PGC = PGC.map(e => mutGen(poi.sample) ::: e) ::: PGC.map(e => mutGen(poi.sample) ::: e)
if (i == pgcRepMale) {
PGCmA = PGC.toArray
}
if (i == pgcRepFemale){
 PGCfA = PGC.toArray
}
}

PGC = Nil

var sgona : List[List[Int]] = Nil
var ogona : List[List[Int]] = Nil

for (i <- 1 to sprmgon) sgona = PGCmA(rand.nextInt(PGCmA.size)) :: sgona
for (i <- 1 to oogon) ogona = PGCfA(rand.nextInt(PGCfA.size)) :: ogona

for (i <- 1 to yold){
for (i <- 1 to 23 ) sgona = sgona.map(e => mutGen(spermpoi.sample) ::: e)
}

println("Meiotic division begins")

var tmpOocyte : List[List[Int]] = Nil
for (x <- ogona){
var a, b: List[Int] = Nil
for (i <- x){	
	if (scala.util.Random.nextInt(2) == 0) a = i :: a //else b = i :: b
}
a = mutGen(1)(0) :: a
tmpOocyte = a :: tmpOocyte
}

println("tmpOocyte " + tmpOocyte.size + " ogona " + ogona.size )

ogona = tmpOocyte

var tmpSperm : List[List[Int]] = Nil

for (x <- sgona){
var a, b: List[Int] = Nil
for (i <- x){	
	if (scala.util.Random.nextInt(2) == 0) a = i :: a else b = i :: b
}
a = mutGen(1)(0) :: a
b = mutGen(1)(0) :: b
tmpSperm = a :: b :: tmpSperm
}
println("tmpSperm " + tmpSperm.size + " sgona " + sgona.size )

sgona = tmpSperm

for (i <- blood){
for (muts <- i){
if (mutFreqBlood.contains(muts)) mutFreqBlood(muts) += 1
else mutFreqBlood += muts -> 1
}
}
val outAFblood = new BufferedWriter(new FileWriter(s"${n}-AF-soma-blood-5y.tab"))
for (d <- mutFreqBlood.toArray.sortBy( a => -1.0 * a._2)){
outAFblood.write(s"${d._2/mutFreqBlood.size.toFloat}\n")
}
outAFblood.close

var AFSperm = new scala.collection.mutable.HashMap[Int, Int]

for (i <- sgona){
for (d <- i){
if (AFSperm.contains(d)) AFSperm(d) += 1
else AFSperm += d -> 1
}
}

val numSperm = sgona.size.toFloat

val outAFsperm = new BufferedWriter(new FileWriter(s"${n}-AF-sperm-5y.tab"))
for (d <- AFSperm.toArray.sortBy( a => -1.0 * a._2)){
outAFsperm.write(s"${d._2/numSperm}\n")
}
outAFsperm.close

val out1Ksperm = new BufferedWriter(new FileWriter(s"${n}-1K-sperm-5y.tab"))
val out1KspermBlood = new BufferedWriter(new FileWriter(s"${n}-1K-sperm-Blood-AF-5y.tab"))

for (i <- 1 to 1000){
for (s <- sgona(rand.nextInt(sgona.size))){
out1Ksperm.write(AFSperm(s)/numSperm + "\t")
if (mutFreqBlood.contains(s)){
out1KspermBlood.write(s"${mutFreqBlood(s)/mutFreqBlood.size.toFloat}\t")
} else {
out1KspermBlood.write("0\t")
}
}
out1Ksperm.write("\n")
out1KspermBlood.write("\n")
}
out1Ksperm.close
out1KspermBlood.close

/* OOcyte */

var AFOocyte = new scala.collection.mutable.HashMap[Int, Int]

for (i <- ogona){
for (d <- i){
if (AFOocyte.contains(d)) AFOocyte(d) += 1
else AFOocyte += d -> 1
}
}

val outAFOocyte = new BufferedWriter(new FileWriter(s"${n}-AF-oocyte-5y.tab"))
val out1KoocyteBlood = new BufferedWriter(new FileWriter(s"${n}-1K-oocyte-Blood-AF-5y.tab"))
for (d <- AFOocyte.toArray.sortBy( a => -1.0 * a._2)){
outAFOocyte.write(s"${d._2/oogon.toFloat}\n")
}
outAFOocyte.close

val out1Koocyte = new BufferedWriter(new FileWriter(s"${n}-1K-oocyte-5y.tab"))
for (i <- 1 to 1000){
for (s <- ogona(rand.nextInt(ogona.size))){
out1Koocyte.write(AFOocyte(s)/oogon.toFloat + "\t")
if (mutFreqBlood.contains(s)){
out1KoocyteBlood.write(s"${mutFreqBlood(s)/mutFreqBlood.size.toFloat}\t")
} else {
out1KoocyteBlood.write("0\t")
}
}
out1KoocyteBlood.write("\n")
out1Koocyte.write("\n")
}
out1Koocyte.close
out1KoocyteBlood.close


}