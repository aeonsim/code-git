//JAVA_OPTS=-Xmx4g  scala -cp ~/Dropbox/PHD-WORK/pedigreeFilter/lib/commons-math3-3.3.jar

object develSims {

import java.io._
import scala.collection.mutable.HashSet

def main (args: Array[String]) : Unit = {
System.err.println("scala app <numKids> <Early Multiplier> <standadRate x.x> <spermRate x.x> <PGCs> <bias?> <earlyDivisions>")

val rand = scala.util.Random
val genomeSize = 2147000000
//val genomeSize = 2600000000L
val addSharedSoma = 3

/* Different differentiation events in embryo development */

val innerCellForm = 7
val innerCellPercent = 0.21
val ebiBlastFormation = 4
val epiBlastPercent = 0.5
val embEbpi = 4
val embEbpiPercent = 0.5

val PGCnum = args(4).toInt
val bias = if (args(5).toUpperCase.contains("T")) true else false

val totalDiv = 36
val lastFastRep = args(6).toInt
val lateDiv = totalDiv - lastFastRep
val mutMulti = args(1).toDouble

val earlyMutRate = (28.0 * ((mutMulti * lastFastRep)/((mutMulti * lastFastRep) + lateDiv)).toDouble)/lastFastRep.toDouble
val lateMutRate = (28.0 * (lateDiv/((args(1).toDouble * lastFastRep) + lateDiv)).toDouble)/lateDiv
System.err.println(s"Early Rate:${earlyMutRate}\nLate Rate:${lateMutRate}")

val offspring = args(0).toInt
val probFilterRange = (1/scala.math.pow(0.5,offspring)).toInt
val pgcRepMale = 18
val pgcRepFemale = 18
val sprmgon = 500000
val oogon = 500000
val yold = 5
val earlypoi = new org.apache.commons.math3.distribution.PoissonDistribution(earlyMutRate)
val poi = new org.apache.commons.math3.distribution.PoissonDistribution(lateMutRate)
val spermpoi = new org.apache.commons.math3.distribution.PoissonDistribution(args(3).toDouble)
var n = 0

while (n < 100){

val mutFreqBlood = new scala.collection.mutable.HashMap[Int, Long]
val mutFreqFemaleGamete = new scala.collection.mutable.HashMap[Int, Int]
val mutFreqMaleGamete = new scala.collection.mutable.HashMap[Int, Int]
val usedMuts = new HashSet[Int]
//for (n <-  (1 to 100)){
 n += 1

def mutGen (num: Int): List[Int] = {
var count = 0
var result : List[Int] = Nil
while (count < num) {
	var mut = rand.nextInt(genomeSize)
	while (usedMuts.contains(mut)) mut = rand.nextInt(genomeSize)
	usedMuts += mut
result = mut :: result
count += 1
}
result
}

//var cells: List[List[Int]] = List(mutGen(poi.sample))
var cells: List[List[Int]] = List(mutGen(0))

/*
if (lastFastRep >= 7) {
// Greater than 7 so first 7 cell divisions high Error rate
for (i <- 1 to 7) cells = cells.map(e => mutGen(earlypoi.sample) ::: e).zip(cells.map(e => mutGen(earlypoi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
} else {
	//not greater so do the high then follow up with the low error rate reps.
for (i <- 1 to lastFastRep) cells = cells.map(e => mutGen(earlypoi.sample) ::: e).zip(cells.map(e => mutGen(earlypoi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
for (i <- lastFastRep to 7) cells = cells.map(e => mutGen(poi.sample) ::: e).zip(cells.map(e => mutGen(poi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
}

var newCells: List[List[Int]] = Nil

//if biased select more closely related cells
for (i <- 1 to (cells.size * innerCellPercent).toInt) {
//if (bias) newCells = cells(i) :: newCells else newCells = cells(rand.nextInt(cells.size)) :: newCells
newCells = cells(rand.nextInt(cells.size)) :: newCells
}

cells = newCells
newCells = Nil

for (i <- 1 to ebiBlastFormation) {
	cells = cells.map(e => mutGen(if (lastFastRep >= 11) earlypoi.sample else poi.sample) ::: e).zip(cells.map(e => mutGen(if (lastFastRep >= 11) earlypoi.sample else poi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
}

for (i <- 1 to (cells.size * epiBlastPercent).toInt) newCells = cells(rand.nextInt(cells.size)) :: newCells

cells = newCells
newCells = Nil

for (i <- 1 to embEbpi) cells = cells.map(e => mutGen(if (lastFastRep >= 15) earlypoi.sample else poi.sample) ::: e).zip(cells.map(e => mutGen(if (lastFastRep >= 15) earlypoi.sample else poi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))

for (i <- 1 to (cells.size * embEbpiPercent).toInt) newCells = cells(rand.nextInt(cells.size)) :: newCells

cells = newCells
newCells = Nil

for (i <- 1 to addSharedSoma) cells = cells.map(e => mutGen(if (lastFastRep >= 18) earlypoi.sample else poi.sample) ::: e).zip(cells.map(e => mutGen(if (lastFastRep >= 18) earlypoi.sample else poi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))

var PGC: List[List[Int]] = Nil

	if (bias){
			for (i <- 1 to PGCnum) PGC = cells(rand.nextInt((cells.size.toDouble * 0.01).toInt)) :: PGC

		} else {
			for (i <- 1 to PGCnum) PGC = cells(rand.nextInt(cells.size)) :: PGC
	}
*/

var crep = 1
var PGC: List[List[Int]] = Nil
var newCells: List[List[Int]] = Nil

while (crep <= 18){

	if (crep <= lastFastRep){

			cells = cells.map(e => mutGen(earlypoi.sample) ::: e).zip(cells.map(e => mutGen(earlypoi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
		
		} else {
		
			cells = cells.map(e => mutGen(poi.sample) ::: e).zip(cells.map(e => mutGen(poi.sample) ::: e)).flatMap(pair => Seq(pair._1,pair._2))
		
		}

	if (crep == 7){
		var grab = 0
		while (grab < (cells.size * innerCellPercent).toInt){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			grab += 1
		}
		cells = newCells
		newCells = Nil
	}

	if (crep == 11){
		var grab = 0
		while (grab < (cells.size * epiBlastPercent).toInt){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			grab += 1
		}
		cells = newCells
		newCells = Nil		
	}

	if (crep == 15){
		var grab = 0
		while (grab < (cells.size * embEbpiPercent).toInt){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			grab += 1
		}
		cells = newCells
		newCells = Nil		
	}

	if (crep == 18){
		println("# of Cells = " + cells.size)
		if (bias){
			for (i <- 1 to PGCnum) PGC = cells(rand.nextInt((cells.size.toDouble * 0.01).toInt)) :: PGC
		} else {
			for (i <- 1 to PGCnum) PGC = cells(rand.nextInt(cells.size)) :: PGC
		}
	}
	crep += 1
}



//for (i <- 1 to 34) PGC = cells(rand.nextInt(cells.size)) :: PGC
val blood = cells
cells = Nil
println(s"After ${addSharedSoma} additional reps have ${PGC.size} PGCs")

var PGCmA, PGCfA = Array(List(0))

var pgcRep = 1
while (pgcRep <= scala.math.max(pgcRepMale,pgcRepFemale)){
//for (i <- 1 to scala.math.max(pgcRepMale,pgcRepFemale)) {
PGC = PGC.map(e => mutGen(poi.sample) ::: e) ::: PGC.map(e => mutGen(poi.sample) ::: e)

if (pgcRep == pgcRepMale) {
PGCmA = PGC.toArray
}
if (pgcRep == pgcRepFemale){
 PGCfA = PGC.toArray
}
pgcRep += 1
}

PGC = Nil

var sgona : List[List[Int]] = Nil
var ogona : List[List[Int]] = Nil


sgona = scala.util.Random.shuffle(PGCmA.toList).take(sprmgon)
ogona = scala.util.Random.shuffle(PGCfA.toList).take(oogon)

var years = 0
while( years < yold){
	var sreps = 0
	while (sreps < 23){
		sgona = sgona.map(e => mutGen(spermpoi.sample) ::: e)
		sreps += 1
	}
	years += 1 
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

val bloodA = blood.toArray
var bcount = 0
println("Blood Cell set Size " + bloodA.size)
//for (i <- blood){
while (bcount < bloodA.size){
	var mutCount = 0	
	//for (muts <- i){
	while (mutCount < bloodA(bcount).size){

		if (mutFreqBlood.contains(bloodA(bcount).apply(mutCount))) {
				mutFreqBlood(bloodA(bcount).apply(mutCount)) += 1L
		} else { 
				mutFreqBlood += bloodA(bcount).apply(mutCount) -> 1L
		}
		mutCount += 1
	}
	bcount += 1
}

//val outAFblood = new BufferedWriter(new FileWriter(s"${n}-AF-soma-blood-5y.tab"))
//for (d <- mutFreqBlood.toArray.sortBy( a => -1.0 * a._2)){
//outAFblood.write(s"${d._1}\t${d._2/mutFreqBlood.size.toFloat}\n")
//}

//outAFblood.close

var AFSperm = new scala.collection.mutable.HashMap[Int, Int]

for (i <- sgona){
for (d <- i){
if (AFSperm.contains(d)) AFSperm(d) += 1
else AFSperm += d -> 1
}
}

val numSperm = sgona.size.toFloat

//val outAFsperm = new BufferedWriter(new FileWriter(s"${n}-AF-sperm-5y.tab"))
//for (d <- AFSperm.toArray.sortBy( a => -1.0 * a._2)){
//outAFsperm.write(s"${d._1}\t${d._2/numSperm}\n")
//}
//outAFsperm.close

val out1Ksperm = new BufferedWriter(new FileWriter(s"${n}-1K-sperm-5y.tab"))
val out1KspermBlood = new BufferedWriter(new FileWriter(s"${n}-1K-sperm-Blood-AF-5y.tab"))

for (i <- 1 to 1000){
	var keep, discard, bloodkeep, blooddiscard : List[String] = Nil
for (s <- sgona(rand.nextInt(sgona.size))){
	if (rand.nextInt(probFilterRange) == 0 ) discard = s"${s}:${AFSperm(s)/numSperm}" :: discard else keep = s"${s}\t${AFSperm(s)/numSperm}" :: keep
	if (mutFreqBlood.contains(s)){
		if (rand.nextInt(probFilterRange) == 0 ) blooddiscard = s"${s}:${mutFreqBlood(s)/mutFreqBlood.size.toFloat}" :: blooddiscard else bloodkeep = s"${s}\t${mutFreqBlood(s)/mutFreqBlood.size.toFloat}" :: bloodkeep
	//out1KspermBlood.write(s"${s}\t${mutFreqBlood(s)/mutFreqBlood.size.toFloat}\t")
	} else {
	//out1KspermBlood.write("${s}\t0\t")
	if (rand.nextInt(probFilterRange) == 0 ) blooddiscard = s"${s}:0" :: blooddiscard else bloodkeep = s"${s}\t0" :: bloodkeep
	}
}
//out1Ksperm.write( s + "\t" + AFSperm(s)/numSperm + "\t")
for (i <- discard.reverse) out1Ksperm.write(i + ",")
if (blooddiscard.size == 0) out1Ksperm.write("-")
out1Ksperm.write("\t")
for (i <- keep.reverse) out1Ksperm.write(i + "\t")

for (i <- blooddiscard.reverse) out1KspermBlood.write(i + ",")
if (blooddiscard.size == 0) out1KspermBlood.write("-")
out1KspermBlood.write("\t")
for (i <- bloodkeep.reverse) out1KspermBlood.write(i + "\t")


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

//val outAFOocyte = new BufferedWriter(new FileWriter(s"${n}-AF-oocyte-5y.tab"))
val out1KoocyteBlood = new BufferedWriter(new FileWriter(s"${n}-1K-oocyte-Blood-AF-5y.tab"))
//for (d <- AFOocyte.toArray.sortBy( a => -1.0 * a._2)){
//outAFOocyte.write(s"${d._1}\t${d._2/oogon.toFloat}\n")
//}
//outAFOocyte.close

val out1Koocyte = new BufferedWriter(new FileWriter(s"${n}-1K-oocyte-5y.tab"))
for (i <- 1 to 1000){
	var keep, discard, bloodkeep, blooddiscard : List[String] = Nil
for (s <- ogona(rand.nextInt(ogona.size))){
	if (rand.nextInt(probFilterRange) == 0 ) discard = s"${s}:${AFOocyte(s)/oogon.toFloat}" :: discard else keep = s"${s}\t${AFOocyte(s)/oogon.toFloat}" :: keep
//out1Koocyte.write(s + "\t" + AFOocyte(s)/oogon.toFloat + "\t")
if (mutFreqBlood.contains(s)){
if (rand.nextInt(probFilterRange) == 0 ) blooddiscard = s"${s}:${(mutFreqBlood(s)/bloodA.size.toDouble)/2.0}" :: blooddiscard else bloodkeep = s"${s}\t${(mutFreqBlood(s)/bloodA.size.toDouble)/2.0}" :: bloodkeep
//out1KoocyteBlood.write(s"${s}\t${mutFreqBlood(s)/mutFreqBlood.size.toFloat}\t")
} else {
if (rand.nextInt(probFilterRange) == 0 ) blooddiscard = s"${s}:0" :: blooddiscard else bloodkeep = s"${s}\t0" :: bloodkeep
}
}
for (i <- discard.reverse) out1Koocyte.write(i + ",")
if (blooddiscard.size == 0) out1Koocyte.write("-")
out1Koocyte.write("\t")
for (i <- keep.reverse) out1Koocyte.write(i + "\t")

for (i <- blooddiscard.reverse) out1KoocyteBlood.write(i + ",")
if (blooddiscard.size == 0) out1KoocyteBlood.write("-")
out1KoocyteBlood.write("\t")
for (i <- bloodkeep.reverse) out1KoocyteBlood.write(i + "\t")




out1KoocyteBlood.write("\n")
out1Koocyte.write("\n")
}
out1Koocyte.close
out1KoocyteBlood.close


}

}
}