import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

val in_ped = new BufferedReader(new FileReader(new File("~/kos_home/refs/CRV-pedigree-800K.ped")))
val in_muts = new BufferedReader(new FileReader(new File("indvs.txt")))


var pedFile = new HashMap[String, Array[String]]
var ofInterest : List[String] = Nil

while (in_muts.ready){
	ofInterest = in_muts.readLine :: ofInterest
}

  /* Recursive function to find all parents & ancestors of nominated individual
	 * returns list of those present in both the PED & VCF files.
	 */

  def itPed(pop: HashMap[String, Array[String]], rec: String, maxIT: Int): List[String] = {
    if (pop.contains(rec) && maxIT >= 0) {
      return rec :: itPed(pop, pop(rec).apply(2), maxIT - 1) ::: itPed(pop, pop(rec).apply(3), maxIT - 1)
    } else {
      return Nil
    }
  }


    var pedigreeList: List[Array[String]] = Nil

    while (in_ped.ready) {
      val temp = in_ped.readLine().split("\t")
      pedFile += temp(1) -> temp
      pedigreeList = temp :: pedigreeList
    }

    in_ped.close


   val ancs = ofInterest.par.map(s => itPed(pedFile,s,40))

   var allIndv = new HashMap[String,Int]

   for (i <- ancs){
   	for (j <- i){
   		if (allIndv.contains(j)) allIndv(j) += 1 else allIndv += j -> 0
   	}
   }

for ( i <- allIndv){
	println( i._1 + "\t" + i._2 + "\t" + i._2/allIndv.size.toFloat )
}
