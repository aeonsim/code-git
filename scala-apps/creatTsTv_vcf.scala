import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

def itPed (pop: HashMap[String, Array[String]] , rec: String, maxIT: Int): List[String] = {
		if (pop.contains(rec) && maxIT >= 0){
			return rec :: itPed(pop,pop(rec).apply(2),maxIT -1) ::: itPed(pop,pop(rec).apply(3),maxIT -1)
		} else {
			return Nil
		}
}

def findChildren (pop: List[Array[String]], vcf: List[String] , parent: String) : List[String]={
		var result: List[String] = Nil
		for (item <- pop){
			if (vcf.contains(item(1)) && ( (item(2) == parent)||(item(3) == parent) )){
				result = item(1) :: result
			}//eif

		}//efor
		result
}//edef


def isHet(genotype:String): Boolean = {
		if (genotype.size == 3) {
			if ( genotype(0) != genotype(2)) true
			else false
		} else false
}

	/* Is this a Variant ie not 0/0 */

def isVar(genotype:String): Boolean = {
		if (genotype.size == 3 && (genotype(0) != '.')) {
			if ((genotype(0) == '0' && genotype(2) == '0')) {
				false
			} else {
 				true
			}
		} else {
 			false
 		}
}
	
	
	/* Is this a Variant ie not ./. */

def isValid(genotype:String): Boolean = {
		if (genotype(0) == '.') {
			false
		} else {
 			true
 		}
}


val in_ped = new BufferedReader(new FileReader(new File(args(1))))
val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(new File(args(0))))))
val out_vcf = new BufferedWriter(new OutputStreamWriter(new BlockCompressedOutputStream(new File("mod_vcf.vcf.gz"))))

var vcfLine = in_vcf.readLine

while (vcfLine(1) == '#'){
	out_vcf.write(vcfLine + "\n")
	vcfLine = in_vcf.readLine
} 

val columnsMap = new HashMap[String, Int]
val columnsMapID = new HashMap[Int, String]
var animals: List[String] = Nil

for (i <- vcfLine.split("\t")){
	val indx = vcfLine.split("\t").indexOf(i)
	columnsMap += i -> indx
	columnsMapID += indx -> i
	if ( indx >= 8 ) animals = i :: animals
}
out_vcf.write(vcfLine + "\n")
/*
* Build Map file of Pedigree
*/

var pedFile = new HashMap[String, Array[String]]

var pedigreeList : List[Array[String]] = Nil

while (in_ped.ready){
		val temp = in_ped.readLine().split("\t")
		pedFile += temp(1) -> temp
		pedigreeList = temp :: pedigreeList
}

	in_ped.close

val families = Array("NL288458773","NL267363937","NL286575328","NL368790870")

var pedData : List[Tuple5[String,List[String],List[String],List[String],HashSet[String]]] = Nil

for (x <- families){
	if (pedFile.contains(x) && animals.contains(x)){
		var parents = List(pedFile(x)(2),pedFile(x)(3))
		var ancestors = itPed(pedFile,x,6).tail.filterNot(x => parents.contains(x) || (! animals.contains(x)))
		var children = findChildren(pedigreeList,animals,x)
		var relatives = new HashSet[String]
		for (y <- animals){
			var ancesL = itPed(pedFile,y,10)

			for (z <- ancestors){
				if (ancesL.contains(z)) relatives += y	
			}
			
		}
		pedData = (x,ancestors,parents,children,relatives) :: pedData
	}
	
}

vcfLine = in_vcf.readLine

while (in_vcf.ready){
	var denovoFake = false
	var hetPro = 0
	var workingPro = ""
	var cLine = vcfLine.split("\t")
	var proband, kids, parents, ances, related, unrelated = 0

	for (p <- families){
		if (isHet(cLine(columnsMap(p)).substring(0,3))) {
			hetPro += 1
			workingPro = p
		} 
	}

	if (hetPro == 1){
		val workingFamily = pedData.filter(_._1 == workingPro)(0)
		//print(pedData.filter(_._1 == workingPro).size)
		for ( i <- 9 to cLine.size -1){
			val cAn = columnsMapID(i)
			//println(workingFamily._1)
			if (isVar(cLine(i).substring(0,3))) {
				cAn match {
					case y if (workingFamily._1 == y) =>  proband += 1
					case y if (workingFamily._3.contains(y)) => parents += 1
					case y if (workingFamily._2.contains(y)) => ances += 1
					case y if (workingFamily._4.contains(y)) => kids += 1
					case y if (workingFamily._5.contains(y)) => related += 1
					case _ => unrelated += 1 
				}

			}

			}


		if (unrelated == 0 && kids >= 1 && parents >= 1){
			for (par <- workingFamily._3){
				//println("candidate")
				print( "cand "+ par + " " + cLine(columnsMap(par)) )
				var tmp = cLine(columnsMap(par)).split(":")
				tmp(0) = "0/0"
				tmp(1) = "25,0"
				tmp(tmp.size - 1) = "0,200,400"
				cLine(columnsMap(par)) = tmp.reduceLeft((a,b) => s"${a}:${b}")
				print(" " + cLine(columnsMap(par)) + "\n")
			}

			for (gp <- workingFamily._2){
				var tmp = cLine(columnsMap(gp)).split(":")
				tmp(0) = "0/0"
				tmp(1) = "25,0"
				tmp(tmp.size - 1) = "0,200,400"
				cLine(columnsMap(gp)) = tmp.reduceLeft((a,b) => s"${a}:${b}")
			}

		}
		}



	//println(s"hetPro:${hetPro} Proband:${proband} Kids:${kids} Unrelated:${unrelated} Anc:${ances} PArents:${parents} Related:${related}")

	if (hetPro == 1 && kids >= 1 && unrelated == 0 && parents >= 1) {
		out_vcf.write(cLine.reduceLeft((a,b) => s"${a}\t${b}") + "\n") 
		println(s"${cLine(0)}\t${cLine(1)}\tModifed Record\t${workingPro}")
		} else {
			out_vcf.write(vcfLine + "\n")
		}

	vcfLine = in_vcf.readLine
}
out_vcf.close
in_vcf.close

