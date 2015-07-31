/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

object analysePhase{

import java.io._
import net.sf.samtools.util._
import scala.collection.mutable.HashMap

/* Global Var for VCF Type*/
var vcfType = ""
val errors = System.err
var PLexist = false
var PL = -1

/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/


	def itPed (pop: HashMap[String, Array[String]] , rec: String): List[String] = {
		if (pop.contains(rec)){
			return rec :: itPed(pop,pop(rec).apply(2)) ::: itPed(pop,pop(rec).apply(3))
		} else {
			return Nil
		}
	}

/* Takes Parental Genotypes & Generates all combinations*/

	def permu(str: String, str1: String): List[String] = {
		val result: List[String] = for (i <- str.toList; j <- str1.toList) yield i.toString.concat(j.toString)
		result
	}

/* Finds Children for a nominated parent based on Pedigree records
* Only finds ones that are both in the VCF & Ped file
*/

	def findChildren (pop: HashMap[String, Array[String]], vcf: HashMap[String,Int] , parent: String) : List[String]={
		var result: List[String] = Nil
		for (item <- pop){
			if (vcf.contains(item._1) && ( (item._2(2) == parent)||(item._2(3) == parent) )){
				result = item._1 :: result
			}//eif

		}//efor
		result
	}//edef

/* Is this a Het variant? */

	def isHet(genotype:String): Boolean = {
		if (genotype.size == 3) {
			if ( genotype(0) != genotype(2)) true
			else false
		} else false
	}

/* Calc Ratio of Ref/Alt*/

	def rc(ref: Int, alt: Int): Float = {
		alt/ref.toFloat
	}

/* Are there two or more Alt reads? */
	def sigAD(alt: Int): Boolean = {
		if (alt >= 2 ) true
		else false
	}

/* Is this a Variant ie not ./. or 0/0 */

	def isVar(genotype:String): Boolean = {
		if (genotype.size == 3 && (genotype(0) != '.')) {
			if ((genotype(0) != '0' && genotype(2) != '0')) {
				true
			} else {
 				false
			}
		} else {
 			false
 		}
	}
	
	/* Return Genotypes as List*/
	
	def rtGTs(genotype: String) : List[String] = {
		val gts = if (genotype.size == 3 && genotype.contains('/')) genotype.split('/') else genotype.split('|')
		gts.toList
	}

	def selROvAD(indv: Array[String], ADval: Int, ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var refAlt = (1,-1)
		vcfType match {
			case "gatk" => refAlt = gatkAD(indv,ADval,GTval)
			case "platypus" => refAlt = platypusRefAlt(indv,ROval, AOval, GTval)
			case "freebayes" => refAlt = fbRefAlt(indv,ROval, AOval, GTval)
			case _ => println("\n\n\nUnknown VCF type in SelROvAD"); System.exit(1)
		}
		refAlt
	}

	def gatkAD(indv: Array[String], ADval: Int, GTval: Int): Tuple2[Int,Int] = {
		if (ADval != -1 && indv.size > ADval && indv(ADval) != "."){
			val RefAlt = indv(ADval).split(",")
			val GT = rtGTs(indv(GTval)).sorted
			(RefAlt(GT(0).toInt).toInt,RefAlt(GT(1).toInt).toInt)
		} else {
		if (PLexist){
			//Check other Genotypes are atleast Phred 10 unlikely
			if (indv(PL).split(",").sorted.tail.forall(_.toInt > 10)) (1,0) else (1,-1)
		}else {
			(1,-1)
			}
		}
	}
	
	def platypusRefAlt(indv: Array[String], ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var alt = -1
		var ref = 1
		
		if (indv.size > AOval){
			val GT = rtGTs(indv(GTval)).sorted
			var refGT = GT(0).toInt
			var altGT = GT(1).toInt
			val alts = indv(AOval).split(",")
			if (ROval != -1 && AOval != -1){
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
				if (refGT == 0){
					ref =  indv(ROval).split(",")(0).toInt - alt
				} else {
					ref = (indv(ROval).split(",")(refGT - 1).toInt - alt) 
				}
			}
		}
		(ref,alt)
	}

	def fbRefAlt(indv: Array[String], ROval: Int, AOval: Int, GTval: Int): Tuple2[Int,Int] = {
		var alt = -1
		var ref = 1
		
		if (indv.size > AOval){
			val GT = rtGTs(indv(GTval)).sorted                                                                                                                                                  
			var refGT = GT(0).toInt
			var altGT = GT(1).toInt
			val alts = indv(AOval).split(",")
			if (ROval != -1 && AOval != -1){
				if (altGT == 0){
					alt = alts(0).toInt
				} else {
					alt = alts(altGT - 1).toInt
				}
				if (refGT == 0){
					ref =  indv(ROval).toInt
				} else {
					ref = alts(refGT - 1).toInt
				}
			}

		} 
		(ref,alt)
	}



/* Take a Genotype String and check against DP limits*/

	def checkDP (genos: Array[String], DPpos: Int, minDP: Int, maxDP: Int): Boolean = {
		if (genos.size >= 2 && DPpos != -1 && genos(DPpos) != "."){
			val curDP = if (vcfType == "platypus") genos(DPpos).split(",")(0).toInt else genos(DPpos).toInt
			if (curDP >= minDP && curDP <= maxDP) true
			else false
		} else {
			false
		}
	}
	
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: Array[String], sire: Array[String], dam: Array[String]) : Tuple2[String, String] = {
	val indvGT = indv(0)
	val sireGT = sire(0)
	val damGT = dam(0)

if (true){
	if(sireGT == "0|0" && damGT == "1|1" && (indvGT == "0|1" || indvGT == "1|0")){
		("0","1")
	} else {
	 if (sireGT == "1|1" && damGT == "0|0" && (indvGT == "0|1" || indvGT == "1|0")){
	 	("1","0")
	 } else {
		if(sireGT == "0|1" && indvGT == "0|1" && damGT == "0|0") {
			("1","0")
		} else {
				if(sireGT == "0|1" && indvGT == "0|1" && damGT == "1|1") {
					("0","1")
				} else {
						if (damGT == "0|1" && indvGT == "0|1" && sireGT == "0|0" ) {
							("0","1")
						} else {
							if (damGT == "0|1" && indvGT == "0|1" && sireGT == "1|1") {
								("1","0")
							} else {
							("x","x")
							}
					}
			}
		}		
	}
	
	}
	}else {("x","x")}

}

/*
def advPhase(curPhase: Tuple2[String,String], child: Array[String], family: Array[String], pedigree: HashMap[String, Array[String]]): String ={
	if (pedigree.contains())
	}
*/

def childPhase(curPhase: Tuple2[String,String], child: Array[String]): String ={
	val childGT = child(0)
	if ((curPhase._1 != "x") && (true) && (childGT == "1|1" || childGT == "0|0")){
		if (curPhase._1 == childGT(0).toString || curPhase._1 == childGT(2).toString) "S" else "D"
		} else {
		"U"
		}
}
	

/* Main body of Program */

def main (args: Array[String]): Unit = {

	var settings = new HashMap[String,String]
	
	for (items <- args){
		val keyVal = items.split("=")
		settings += keyVal(0).toUpperCase -> keyVal(1) 
	}
	

	if ((! settings.contains("VCF")) && (! settings.contains("PED")) & (! settings.contains("TRIOS"))) {
		println("advFilter VCF=input.vcf.gz PED=input.ped TRIOS=input_probands.txt ")
		println("Trios = txtfile per line: AnimalID\tavgDepth")
		println("{} Optional arguments, Values shown are default")
		System.exit(1)
	}

	val in_vcf = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream(settings("VCF")))))
	val in_ped = new BufferedReader(new FileReader(settings("PED")))
	val in_pro = new BufferedReader(new FileReader(settings("TRIOS")))
		
	/* Default settings*/
	
	var pedFile = new HashMap[String, Array[String]]

	// Tuple6(Ancestors, Parents, Children, Descendants, TrioDP, population, extended-Fam)

	var trios = new HashMap[String, Tuple7[List[String],List[String],List[String],List[String],Int,List[String],List[String]]]
	var vcfanimals = new HashMap[String, Int]
	var ancestors : List[String] = Nil
	var parents : List[String] = Nil
	var children : List[String] = Nil
	var descendents : List[String] = Nil
	var extFam : List[String] = Nil
	var population : List[String] = Nil
	var lastPhase = new HashMap[String,Tuple2[String,String]]
	var childState = new HashMap[String,String]
	var allChildren = new HashMap[String,String]
	var childHaps = new HashMap[String, List[String]]
/*
* Iterate through VCF file to find Column headers and store positions in Map for later recall
*/

	var vcfrec = in_vcf.readLine().split("\t")
	while (vcfrec(0).apply(1) == '#'){
		vcfrec = in_vcf.readLine().split("\t")
	}

	if (vcfrec(0).apply(0) == '#' && vcfrec(0).apply(1) == 'C'){
		for (i <- 9 to (vcfrec.size -1)){
			vcfanimals += vcfrec(i) -> i
		}
	}

	val animalIDS = vcfanimals.keys.toArray

/*
* Build Map file of Pedigree
*/

	while (in_ped.ready){
		val temp = in_ped.readLine().split("\t")
		pedFile += temp(1) -> temp
	}

	in_ped.close

/*
* Iterate through our list of probands finding their sequenced Ancestors, Parents, Children
* & other descendants, to use for filtering, and build the Trios Map
*/

	while (in_pro.ready){
		val curPro = in_pro.readLine.split("\t")
		descendents = Nil
		extFam = Nil
		if (pedFile.contains(curPro(0)) && animalIDS.contains(curPro(0))){
			parents = pedFile(curPro(0))(2) :: pedFile(curPro(0))(3) :: Nil
			ancestors = itPed(pedFile,curPro(0)).tail.filterNot(x => parents.contains(x) || (! vcfanimals.contains(x)))
			children = findChildren(pedFile,vcfanimals,curPro(0))
			for (indv <- animalIDS.toList){ //.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)))){
				val tempPed = itPed(pedFile,indv)
				if (tempPed.contains(curPro(0))){
					descendents = indv :: descendents
				}
				if (pedFile.contains(indv) && ( pedFile(indv)(2) == parents(0) || pedFile(indv)(3) == parents(0) )  && pedFile(indv)(1) != curPro(0)){
					extFam = indv :: extFam
				}
				if (! pedFile.contains(indv)){
					pedFile += indv -> Array("noFAM",indv,"0","0","0")
				}
			}//efor
			var tmpdesc = descendents.filterNot(x => ((x == curPro(0)) || children.contains(x) || ancestors.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			population = animalIDS.toList.filterNot(x => ( (x == curPro(0)) || children.contains(x) || ancestors.contains(x) || descendents.contains(x) || parents.contains(x)|| (! vcfanimals.contains(x))))
			
			if (animalIDS.contains(parents(0)) && animalIDS.contains(parents(1)) && animalIDS.contains(curPro(0)) && children.size != 0){
				trios += curPro(0) -> (ancestors, parents, children, tmpdesc, curPro(1).toInt, population, extFam)
				for (child <- children) {
				allChildren += child -> ""
				childHaps += child ->  ("Chr\tStart\tEnd\tState" :: Nil)
				}
			}
		}//eif
	}//ewhile

	in_pro.close

println("Built Pedigrees\n")
//Closed proband file replace println with statsOut.write

/* Report identified Trios & there Pedigree Structure*/

	for (fam <- trios){
		print(s"TRIO:\t${fam._1}\t${fam._2._2(0)}\t${fam._2._2(1)}\n")
		print("Grandparents\t")
		fam._2._1.foreach(s => println(s + "\t"))
		print(fam._2._1.size + "\n")
		print("Children\t")		
		fam._2._3.foreach(s => println(s + "\t"))
		print(fam._2._3.size + "\n")
		print("Descendents\t")
		fam._2._4.foreach(s => println(s + "\t"))
		print(fam._2._4.size + "\n")
		print("EXFAM\t")
		fam._2._7.foreach(s => println(s + "\t"))
		print(fam._2._7.size + "\n\n")
		lastPhase += fam._1 -> ("u","u")
		childState += fam._1 -> ""
	}


	while (in_vcf.ready){
	val curTrio = "NL288458773"
	val cline = in_vcf.readLine.split("\t")
	val fam = trios(curTrio)
	val curSNP = phase(cline(vcfanimals(curTrio)).split(":"),cline(vcfanimals(fam._2(0))).split(":"),cline(vcfanimals(fam._2(1))).split(":"))
	if (curSNP != ("x","x")){
		//fam._2 Parents, fam._3 Children
	//println(fam._2 + "\t" + fam._3)
	print(s"${fam._2(0)}:${cline(vcfanimals(fam._2(0)))}\t${fam._2(1)}:${cline(vcfanimals(fam._2(1)))}\t${curTrio}:${cline(vcfanimals(curTrio))}" )
	print("\t" + curSNP)
	}
	for (child <- fam._3){
	print ("\t" + child + ": ")
	val pars = pedFile(child)
	if (vcfanimals.contains(pars(2)) && vcfanimals.contains(pars(3))){
	val phaseGT = phase(cline(vcfanimals(child)).split(":"),cline(vcfanimals(pars(2))).split(":"),cline(vcfanimals(pars(3))).split(":"))
	if (phaseGT != ("x","x")) print(phaseGT) else print (cline(vcfanimals(child)) + " Shp")
	print("\t")
	} else {
	val cPhase = childPhase(curSNP,cline(vcfanimals(child)).split(":"))
	 print (cPhase + " " + cline(vcfanimals(child)) + " Shp\t")
	}
	
	}
	print("\n")
	
	}

	
}//eMain

}//eObject