import scala.collection.mutable.HashMap
import scala.collection.parallel._

if (args.size == 0) System.err.println("app.scala <baseRate x.x> <early *X> <popsize> <male2female 1.0> <generations> <mosaic true/false>")
if(args.size == 0) System.exit(0)
val genome = 2000000000
val rand = new scala.util.Random

//val baseRate = 0.4
val multiplier = args(1).toInt
val popsize = args(2).toInt
val replications = 14
val male2female = args(3).toDouble
val generations = args(4).toInt
val mosaic = args(5).toUpperCase.contains("T")
var toRemove : List[Int] = Nil

def mutGen (num: Int): List[Int] = {
	val mos = if (num < 0) true else false
	var count = 0
	var result : List[Int] = Nil
		while (count < scala.math.abs(num)) {
			val tmp = if (mos) -1 * rand.nextInt(genome) else rand.nextInt(genome)
			result = tmp :: result
			count += 1
		}
	result
}

//val baseRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble)
//val earlyRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble * multiplier)
//val spermRate = new org.apache.commons.math3.distribution.PoissonDistribution(if (mosaic) 6 else 12)
//Arrays of the Mutations carried by each individual
var curPop = Array.fill(popsize)(List((-1,-1)))
var nextGen = Array.fill(popsize)(List((-1,-1)))
var gametes = new Array[List[List[Int]]] (popsize)
var curGen = 0
val mutAge = new HashMap[Int,Int]
var fixtime : List[Int] = Nil
var mosDegree : Map[Int,Double] = Map()
var mosMutFix : List[Int] = Nil
var nmMutFix : List[Int] = Nil
var avgNumMuts: List[Int] = Nil
var numMos, numNM = 0
var mLast500, nmLast500 : List[Int] = Nil

def gam(l:List[Tuple2[Int,Int]]): List[Tuple2[Int,Int]] ={
      l.map(x => if(x._2 == 2) (x._1,1) else if (rand.nextInt(2) == 1 ) (x._1,1) else (x._1,0)).filter( x => x._2 != 0)

  }

while (curGen < generations){
	var indv = 0
	for (indv2 <- (0 to popsize -1).par){
		var cells: List[List[Int]] = List(mutGen(0))
		var reps = 0
		while (reps < replications) {
			
			if (mosaic){
				if (reps == 1){
					cells = cells.map(e => mutGen(-1) ::: e) ::: cells.map(e => mutGen(-1) ::: e)
					//System.err.println(cells.size)	
				} else {
					//if (reps == (replications - 1) && indv2 < (male2female * popsize)) cells = cells.map(e => mutGen(1) ::: e) ::: cells.map(e => mutGen(1) ::: e) else cells = cells.map(e => mutGen(0) ::: e) ::: cells.map(e => mutGen(0) ::: e)
					// Want only 1 mosaic and non-mosaic mut
					if (reps == (replications - 1) ) cells = cells.map(e => mutGen(1) ::: e) ::: cells.map(e => mutGen(1) ::: e) else cells = cells.map(e => mutGen(0) ::: e) ::: cells.map(e => mutGen(0) ::: e)
				}
			} else {
				if (reps == (replications - 1)){
					if (indv2 < (male2female * popsize)) cells = cells.map(e => mutGen(1) ::: e) ::: cells.map(e => mutGen(1) ::: e) else cells = cells.map(e => mutGen(1) ::: e) ::: cells.map(e => mutGen(1) ::: e)
				} else {
					cells = cells.map(e => mutGen(0) ::: e) ::: cells.map(e => mutGen(0) ::: e)
				}
			}
			reps += 1

		}

		gametes(indv2) = cells
		//indv += 1
	}
	indv = 0
	//If zero Gen muts = gametes, randomly take a gamete from two individuals and merge
	if (curGen == 0){
		//println( curGen + " 1st Gen Gametes")
			//for (indv2 <- (0 to popsize -1).par){
			while(indv < popsize){
				val father = rand.nextInt((male2female * popsize).toInt)
				val mother = rand.nextInt(((1 - male2female ) * popsize).toInt) + (male2female * popsize).toInt
				val maleGam = gametes(father)(rand.nextInt(gametes(father).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1))
				val femaleGam = gametes(mother)(rand.nextInt(gametes(mother).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1))
				nextGen(indv) = (maleGam.toList ::: femaleGam.toList).groupBy(_._1).mapValues(_.map(_._2).sum).toList
				nextGen(indv).foreach(x => if(x._2 > 1)System.err.println(x + "M"))
				indv += 1
			}
		} else {
			while(indv < popsize){
				val father = rand.nextInt((male2female * popsize).toInt)
				val mother = rand.nextInt(((1 - male2female ) * popsize).toInt) + (male2female * popsize).toInt
				val maleGam = (gametes(father)(rand.nextInt(gametes(father).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(father))).filter(x => ! toRemove.contains(x._1)).distinct
				val femaleGam = (gametes(mother)(rand.nextInt(gametes(mother).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(mother))).filter(x => ! toRemove.contains(x._1)).distinct
				nextGen(indv) = (maleGam ::: femaleGam).groupBy(_._1).mapValues(_.map(_._2).sum).toList
				indv += 1
			}
		}
		toRemove = Nil
		//build list of all muts and number of times seen
		val muts = new HashMap[Int,Int]
		indv = 0
		while(indv < popsize){
			var mutnum = 0
			while (mutnum < nextGen(indv).size){
				if (muts.contains(nextGen(indv)(mutnum)._1)){
						if (nextGen(indv)(mutnum)._2 == 2) muts(nextGen(indv)(mutnum)._1) += 1
					} else {
						muts +=  nextGen(indv)(mutnum)._1 -> 1
					}
				mutnum += 1
			}
			indv += 1
		}
		//Check to see if mut is new if not add one gen onto age if is add and set to 1
		for (m <- muts){
			if (mutAge.contains(m._1)){
				mutAge(m._1) += 1
			} else {
				mutAge += m._1 -> 1
			}
			if (m._1 < 0) numMos += 1 else numNM += 1

			if (m._2 == (popsize)) {
				if (m._1 < 0) {
					mosMutFix = mutAge(m._1) :: mosMutFix
					if (mLast500.size == 500){
							mLast500 = mutAge(m._1) :: mLast500.dropRight(1)
						} else {
							mLast500 = mutAge(m._1) :: mLast500
					}
				} else{
					nmMutFix = mutAge(m._1) :: nmMutFix
					if (nmLast500.size == 500){
							nmLast500 = mutAge(m._1) :: nmLast500.dropRight(1)
						} else {
							nmLast500 = mutAge(m._1) :: nmLast500
					}
				} 
				//println(s"${curGen}\t${if( m._1 < 0 )"m" else "nm"}\t${mutAge(m._1)}\t${if(m._1 < 0 ) mosMutFix.size else nmMutFix.size}\t${if(m._1 < 0 ) mosMutFix.sum/mosMutFix.size.toDouble else nmMutFix.sum/nmMutFix.size.toDouble}\t${if(m._1 < 0 ) mosMutFix.sorted.apply(mosMutFix.size/2) else nmMutFix.sorted.apply(nmMutFix.size/2)}" + s"\t${mosMutFix.size} / " + numMos + s"\t${nmMutFix.size} / " + numNM + 
				//	"\t" + mLast500.sum/mLast500.size.toDouble + "\t" + nmLast500.sum/nmLast500.size.toDouble + "\t" + mLast500.sorted.apply(mLast500.size/2) + "\t" + nmLast500.sorted.apply(nmLast500.size/2))
				println(s"${curGen}\t${if( m._1 < 0 )"m" else "nm"}\t${mutAge(m._1)}\t${if(m._1 < 0 ) mosMutFix.size else nmMutFix.size}\t${if(m._1 < 0 ) mosMutFix.sum/mosMutFix.size.toDouble else nmMutFix.sum/nmMutFix.size.toDouble}\t${if(m._1 < 0 ) mosMutFix.sorted.apply(mosMutFix.size/2) else nmMutFix.sorted.apply(nmMutFix.size/2)}" + s"\t${mosMutFix.size} / " + numMos + s"\t${nmMutFix.size} / " + numNM + 
					"\t" + mLast500.sum/mLast500.size.toDouble + "\t" + nmLast500.sum/nmLast500.size.toDouble)

				fixtime = mutAge(m._1) :: fixtime
				toRemove = m._1 :: toRemove
			} else {
			}

		}

		//check to see if mut has reach fixation if so report age and add to ignore

		curPop = nextGen

	curGen += 1
}

System.err.println(mosMutFix.sum/mosMutFix.size.toDouble + "\t" + nmMutFix.sum/nmMutFix.size.toDouble + s"\t${mosMutFix.size} / " + numMos + s"\t${nmMutFix.size} / " + numNM)


