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

val baseRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble)
val earlyRate = new org.apache.commons.math3.distribution.PoissonDistribution(args(0).toDouble * multiplier)
val spermRate = new org.apache.commons.math3.distribution.PoissonDistribution(if (mosaic) 6 else 12)
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
//({ case ((k,v1),(_,v2)) => (k,v1+v2) })

def gam(l:List[Tuple2[Int,Int]]): List[Tuple2[Int,Int]] ={
	  //l.foreach(x => if(x._2 > 2)System.err.println(x))
      l.map(x => if(x._2 == 2) (x._1,1) else if (rand.nextInt(2) == 1 ) (x._1,1) else (x._1,0)).filter( x => x._2 != 0)

  }

 //gam(sperm).merged(gam(egg))({ case ((k,v1),(_,v2)) => (k,v1+v2) })

while (curGen < generations){
	var indv = 0
	//println("Gen Gametes")
	for (indv2 <- (0 to popsize -1).par){
	//while(indv < popsize){
		var cells: List[List[Int]] = List(mutGen(0))
		var reps = 0
		while (reps < replications) {
			
			if (mosaic){
				//println("mos")
				if (reps < 4){
					cells = cells.map(e => mutGen(-1 * earlyRate.sample) ::: e) ::: cells.map(e => mutGen(-1 * earlyRate.sample) ::: e)	
				} else {
					if (reps == (replications - 1) && indv2 < (male2female * popsize)) cells = cells.map(e => mutGen(spermRate.sample) ::: e) ::: cells.map(e => mutGen(spermRate.sample) ::: e) else cells = cells.map(e => mutGen(0) ::: e) ::: cells.map(e => mutGen(0) ::: e)
				}
				//reps += 1
			} else {
				if (reps == (replications - 1)){
					if (indv2 < (male2female * popsize)) cells = cells.map(e => mutGen(spermRate.sample) ::: e) ::: cells.map(e => mutGen(spermRate.sample) ::: e) else cells = cells.map(e => mutGen(baseRate.sample) ::: e) ::: cells.map(e => mutGen(baseRate.sample) ::: e)
				} else {
					cells = cells.map(e => mutGen(baseRate.sample) ::: e) ::: cells.map(e => mutGen(baseRate.sample) ::: e)
				}
				//reps += 1

			}
			reps += 1

		}

		//val mosCount = cells.flatten.groupBy(identity).map(x => (x._1, x._2.size)).map(x => (x._1,x._2/cells.size.toDouble))
		//println("avg muts " + cells.map(_.size).sum/cells.size)
		//val tmp : scala.collection.immutable.HashMap[Int,Double] = mosCount.map(x => (x._1,x._2/cells.size.toDouble))
		//val tmpMD = collection.mutable.HashMap() ++ mosCount ++ mosDegree
		//mosDegree = tmpMD
		//System.err.println(cells.size)
		//cells(rand.nextInt(cells.size)).foreach(x => System.err.println(x))
		//System.err.println(cells(rand.nextInt(cells.size)).size)
		gametes(indv2) = cells
		//indv += 1
	}
	//println("calc ARs")
	//gametes.flatten.flatten.groupBy(identity).map(x => (x._1, x._2.size)).map(x => (x._1,x._2/cells.size.toDouble))
	//mosDegree = gametes.par.map(y => y.flatten.groupBy(identity).map(x => (x._1, x._2.size)).map(x => (x._1,x._2/y.size.toDouble))).flatten.toMap ++ mosDegree
	//println("done calc ARs")
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
			//combines two tuple lists and combines the values list((2,1),(34,2))
			//(z ::: l).groupBy(_._1).mapValues(_.map(_._2).sum)
			//println( curGen + " Gen Gametes")
			//for (indv2 <- (0 to popsize -1).par){
			while(indv < popsize){
				val father = rand.nextInt((male2female * popsize).toInt)
				val mother = rand.nextInt(((1 - male2female ) * popsize).toInt) + (male2female * popsize).toInt
				//val maleGam = (gametes(father)(rand.nextInt(gametes(father).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(father))).filter(x => ! toRemove.contains(x._1)).groupBy(_._1).mapValues(_.map(_._2).sum).toList.map(z => (z._1,1))
				val maleGam = (gametes(father)(rand.nextInt(gametes(father).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(father))).filter(x => ! toRemove.contains(x._1)).distinct
				//maleGam.foreach(x => if(x._2 > 2)System.err.println(x + "F"))
				//val femaleGam = (gametes(mother)(rand.nextInt(gametes(mother).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(mother))).filter(x => ! toRemove.contains(x._1)).groupBy(_._1).mapValues(_.map(_._2).sum).toList.map(z => (z._1,1))
				val femaleGam = (gametes(mother)(rand.nextInt(gametes(mother).size)).filter(y => rand.nextInt(2) == 1).map(y => (y,1)) ::: gam(curPop(mother))).filter(x => ! toRemove.contains(x._1)).distinct
				//femaleGam.foreach(x => if(x._2 > 2)System.err.println(x + "M"))
				nextGen(indv) = (maleGam ::: femaleGam).groupBy(_._1).mapValues(_.map(_._2).sum).toList
				//nextGen(indv).foreach(x => if(x._2 > 2)System.err.println(x + "M"))
				indv += 1
			}
		}
		toRemove = Nil
		//build list of all muts and number of times seen
		val muts = new HashMap[Int,Int]
		indv = 0
		//for (indv2 <- (0 to popsize -1).par){
		while(indv < popsize){
			//println("Counting Muts")
			var mutnum = 0
			while (mutnum < nextGen(indv).size){
				if (muts.contains(nextGen(indv)(mutnum)._1)){
						if (nextGen(indv)(mutnum)._2 == 2) muts(nextGen(indv)(mutnum)._1) += 1
					} else {
						//if (nextGen(indv)(mutnum)._2 == 2) muts +=  nextGen(indv)(mutnum)._1 -> 1
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
				//println(s"Mutation ${m._1} fixed ${m._2} in ${mutAge(m._1)} was at AR ${mosDegree(m._1)}")
				//if (mosDegree(m._1) >= 0.01) mosMutFix = mutAge(m._1) :: mosMutFix else nmMutFix = mutAge(m._1) :: nmMutFix
				if (m._1 < 0) {
					mosMutFix = mutAge(m._1) :: mosMutFix 
				} else{
					nmMutFix = mutAge(m._1) :: nmMutFix
				} 
				println(s"${curGen}\t${if( m._1 < 0 )"m" else "nm"}\t${mutAge(m._1)}\t${if(m._1 < 0 ) mosMutFix.size else nmMutFix.size}\t${if(m._1 < 0 ) mosMutFix.sum/mosMutFix.size.toDouble else nmMutFix.sum/nmMutFix.size.toDouble}\t${if(m._1 < 0 ) mosMutFix.sorted.apply(mosMutFix.size/2) else nmMutFix.sorted.apply(nmMutFix.size/2)}" + s"\t${mosMutFix.size} / " + numMos + s"\t${nmMutFix.size} / " + numNM)
				fixtime = mutAge(m._1) :: fixtime
				toRemove = m._1 :: toRemove
				//mosDegree.remove(m._1)
			} else {
				//if ((m._2 / popsize.toDouble) >= 0.9) System.err.println(s"Mutation ${m._1} at ${m._2/popsize.toDouble} in ${mutAge(m._1)}")
			}

		}

		//check to see if mut has reach fixation if so report age and add to ignore

		curPop = nextGen

	curGen += 1
}

//System.err.println("Average fix time = " + fixtime.sum/fixtime.size.toDouble + " Num of fixed = " + fixtime.size)
//System.err.println("M >= 1% avg fix  = " + mosMutFix.sum/mosMutFix.size.toDouble + " Num of fixed = " + mosMutFix.size.toDouble)
//System.err.println("Median Mos fix time = " + mosMutFix.sorted.apply(mosMutFix.size/2))
//System.err.println("NM <= 1% avg fix = " + nmMutFix.sum/nmMutFix.size.toDouble + " Num of fixed = " + nmMutFix.size.toDouble) 
//System.err.println("Median NM fix time = " + nmMutFix.sorted.apply(nmMutFix.size/2))
//System.err.println("avg number of muts = " + avgNumMuts.sum/avgNumMuts.size.toDouble )
System.err.println(mosMutFix.sum/mosMutFix.size.toDouble + "\t" + nmMutFix.sum/nmMutFix.size.toDouble + s"\t${mosMutFix.size} / " + numMos + s"\t${nmMutFix.size} / " + numNM)
//System.err.println(mosMutFix.sorted.apply(0) + "\t" + nmMutFix.sorted.apply(0))
//System.err.println(curPop.map(_.size).sum/curPop.size)

