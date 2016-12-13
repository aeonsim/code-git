//JAVA_OPTS=-Xmx4g  scala -cp ~/Dropbox/PHD-WORK/pedigreeFilter/lib/commons-math3-3.3.jar

import java.io._

val rand = scala.util.Random
val genomeSize = 2147000000
//val genomeSize = 2600000000L
val addSharedSoma = 3

/* Different differentiation events in embryo development */

val innerCellForm = 7
val innerCellPercent = 0.20
val ebiBlastFormation = 4
val epiBlastPercent = 0.5
val embEbpi = 4
val embEbpiPercent = 0.5


val pgcRepMale = 18
val pgcRepFemale = 18
val sprmgon = 500000
val oogon = 500000
val yold = 5

def mutGen (num: Int): List[Int] = {

var count = 0
var result : List[Int] = Nil
while (count < num) {
val chr = if (rand.nextInt(2) == 0) -1 else 1
result = (chr * rand.nextInt(genomeSize)) :: result
count += 1
}
result
}

System.err.println("scala app <earlymultiplier> <numMuts> <numReps> <bias_t/f> <numPGCs> <numEarlyReps>")

val mult = args(0).toInt
val numMutations = args(1).toInt
val cellCycles = args(2).toInt
val biasedBool = if (args(3).toUpperCase.contains("T")) true else false
val pgcCount = args(4).toInt
val numEarly = args(5).toInt

val pgcSelection = if (pgcCount < 10) 14 else 17

//Randomly assign mutations to replications

var mutations = Array.fill[Int](cellCycles){0}
var PGC: List[List[Int]] = Nil


for (i <- 1 to numMutations) {
mutations(rand.nextInt(cellCycles)) += 1
}


if (mult != 1){
	println("4 rep multiple is ")
	//val multiple = args(0).toInt
	
	mutations = Array.fill[Int](cellCycles){0}

	for (i <- 1 to numMutations) {
		val rnum = cellCycles + ((mult - 1 ) * numEarly)
		val nexMut = rand.nextInt(rnum)
	
		if (nexMut < cellCycles){
	
			mutations(nexMut) += 1
		} else {
	
			val tmpMu = nexMut - cellCycles
	
			val newVal = 1 + (tmpMu/(mult - 1))
	
			mutations(newVal) += 1
		}
		
	}
}



//create replication method

def mutQ(rep: Int, masterMut: Array[Int]): List[Int] ={
	mutGen(masterMut(rep))
}

// Record total reps

var divisions = 0
var masterCell : List[Int] = Nil
var cells, somaCells: List[List[Int]] = Nil


// Begin cell cycle

while (divisions < cellCycles){
	//if (divisions < 40 ) cells = masterCell :: cells ::: cells
	if (divisions < 40 ) cells = masterCell :: cells.zip(cells).flatMap(pair => Seq(pair._1,pair._2))
	masterCell = masterCell ::: mutQ(divisions,mutations)

	if (divisions == 7){
		val keepCells = innerCellPercent * (cells.size + 1)
		var countC = 0
		var newCells: List[List[Int]] = Nil 
		while (countC <= keepCells){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			countC += 1
		}
		println("Inner cell mass formation keeping " + newCells.size + " of " + (cells.size + 1))
		cells = newCells

	}

	if (divisions == 11){
		val keepCells = epiBlastPercent * (cells.size + 1)
		var countC = 0
		var newCells: List[List[Int]] = Nil 
		while (countC <= keepCells){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			countC += 1
		}
		println("Epiblast formation keeping " + newCells.size + " of " + (cells.size + 1))
		cells = newCells
	}

	if (divisions == 15){
		val keepCells = embEbpiPercent * (cells.size + 1)
		var countC = 0
		var newCells: List[List[Int]] = Nil 
		while (countC <= keepCells){
			newCells = cells(rand.nextInt(cells.size)) :: newCells
			countC += 1
		}
		println("Embryonic Epiblast formation keeping " + newCells.size + " of " + (cells.size + 1))
		cells = newCells
	}

	//if (divisions >= 14 && divisions <= 17){
	if (divisions == pgcSelection ){

		if (biasedBool){
			/*var biasArray :List[List[Int]] = Nil

			var cCell  = 0
			while (cCell < cells.size){
				//Add multiplier if what stronger bias
				var nMut = 1 + (cells(cCell).size * cells(cCell).size)
				var muC = 0
				while (muC < nMut){
					biasArray = cells(cCell) :: biasArray
					muC += 1
				}

				cCell += 1
			}


			for (i <- 1 to (pgcCount - 1)) PGC = biasArray(rand.nextInt(biasArray.size)) :: PGC
			*/
			for (i <- 1 to (pgcCount - 1)) PGC = cells(rand.nextInt((cells.size * 0.2).toInt)) :: PGC
		} else {
			for (i <- 1 to (pgcCount - 1)) PGC = cells(rand.nextInt(cells.size)) :: PGC
		}

		}
	if (divisions == 17) {
		println(PGC.size + " sampled")
		println("Primary lineage  at formation of PGCs contains " + masterCell.size + " de novos,  " + masterCell.count(a => a > 0) + " Paternal and " + masterCell.count(a => a < 0) + " Maternal")
		somaCells = cells
		cells = PGC.tail
	}


	if (divisions == 35){
		println("Maximum PGCs reached " + cells.size)
		val keepCells = 1000000
		var newCells: List[List[Int]] = Nil 
		newCells = scala.util.Random.shuffle(cells).take(keepCells)
		println("Downsampled PGCs " + newCells.size + " of " + (cells.size + 1))
		cells = newCells

	}

	//if(divisions > 36 ) print( divisions + " ")


divisions += 1
}

println("Have " + masterCell.size + " de novos after " + divisions + " cell divisions")

var finalCells = cells.toArray
//cells = Nil

//var mutations = Array.fill[Int](cellCycles){0}

var MutHaps = Array.fill[Array[Int]](numMutations + 1){Array.fill[Int](numMutations + 1){0}}

var resCount = 0

while(resCount < finalCells.size){
	var mutCounter = 0
	var mutA, mutB, hapA, hapB = 0
	while (mutCounter < masterCell.size){
		//println(finalCells(resCount))

	if (finalCells(resCount).contains(masterCell(mutCounter))){
			if (rand.nextInt(2) == 0) {
				mutA += 1
				hapA += 1
			} else{
				mutB += 1
				hapB += 1
			}
		} else{
			if (rand.nextInt(2) == 0)  hapA += 1 else hapB += 1
		}

		mutCounter += 1

	}

	MutHaps(mutA)(hapA) += 1
	MutHaps(mutB)(hapB) += 1

	resCount += 1
}

for (h <- 0 to numMutations){
	print("\tHaplo" + h)
}

print("\n")
for (m <- 0 to (numMutations )){
print("Muts_" + m )
for (h <- 0 to (numMutations)){
	print("\t" + MutHaps(m)(h))
}
print("\n")
}
