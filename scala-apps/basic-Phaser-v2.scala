
/* Advanced Pedigree Filter
* Takes in VCF (gz) & Pedigree file, with cutoffs for DP, and a list of Proband ID's
* Use pedigree file to build family tree for each Proband, Proband, Children, Descendants, Parents, Ancestors
* Select variants based on pedigree & min DP check AD for low coverage
* Parents & Ancestors -ve, Proband +ve, Children 1+ positive, Descendants possibly +ve
* Output VCF per Proband ID of possible De novos along with, score file with counts in the population.
* Check to see if each VCF line is the correct length & discard if not.
*/

/* phaseTracker takes in phase blocks and tracks position within the block for each query*/
import scala.collection.mutable.HashMap
import scala.math._

class phaseTracker (phaseData: HashMap[String,Array[Tuple5[String,Int,Int,String,Int]]]){



	
	var curPhaseBlocks = new HashMap[String,Int]
	for (c <- phaseData.keys){
		curPhaseBlocks += c -> 0
	}
	
	
	def getPhase(chrom: String, position: Int) : String ={
	
	//Simple block ID
		val block = phaseData(chrom).filter(pos => (pos._2 <= position && pos._3 >= position))
		if (block.size == 1){
			block(0)._4
		} else {
			"BOTH"
		}

			
	} //ed def
	
	def getNearestBlock(chrom: String, position: Int) : String = {
			val next = phaseData(chrom).filter(pos => (pos._2 >= position))
			val prev = phaseData(chrom).filter(pos => (pos._3 <= position))
			if (next.size >= 1 && prev.size >= 1){
			  val nextDistance = scala.math.abs(position - next.head._2)
			  val prevDistance = scala.math.abs(position - prev.head._3)
			  if (nextDistance > prevDistance) next.head._4 else prev.head._4
			} else {
			  if ((next.size >= 1 && prev.size == 0)) next.head._4
			  else if (next.size == 0 && prev.size >= 1) prev.head._4 else "BOTH"
			}
			
	
	}
}


object simplePhase {

import java.io._
import scala.collection.mutable.HashMap


var GQ, DP = 0

/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/


	def itPed (pop: HashMap[String, Array[String]] , rec: String, maxIT: Int): List[String] = {
		if (pop.contains(rec) && maxIT >= 0){
			return rec :: itPed(pop,pop(rec).apply(2),maxIT -1) ::: itPed(pop,pop(rec).apply(3),maxIT -1)
		} else {
			return Nil
		}
	}

/* Finds Children for a nominated parent based on Pedigree records
* Only finds ones that are both in the VCF & Ped file
*/

	def findChildren (pop: List[Array[String]], vcf: Array[String] , parent: String) : List[String]={
		var result: List[String] = Nil
		for (item <- pop){
			if (vcf.contains(item(1)) && ( (item(2) == parent)||(item(3) == parent) )){
				result = item(1) :: result
			}//eif

		}//efor
		result
	}//edef
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: String, sire: String, dam: String) : Tuple2[String, String] = {
if (sire != "." && dam != "." && indv != "."){
	if(sire == "0/0" && dam == "1/1" && (indv == "1/0" || indv == "0/1")){
		("0","1")
	} else {
	 	if (sire == "1/1" && dam == "0/0" && (indv == "0/1" || indv == "1/0")){
	 		("1","0")
	 		} else {
				if((sire == "0/1" || sire == "1/0") && (indv == "0/1" || indv == "1/0") && dam == "0/0") {
					("1","0")
				} else {
					if((sire == "0/1" || sire == "1/0") && (indv == "0/1" || indv == "1/0") && dam == "1/1") {
						("0","1")
					} else {
						if ((dam == "0/1" || dam == "1/0") && (indv == "0/1" || indv == "1/0") && sire == "0/0" ) {
							("0","1")
						} else {
							if ((dam == "0/1" || dam == "1/0") && (indv == "0/1" || indv == "1/0") && sire == "1/1") {
								("1","0")
							} else {
								("x","x")
								}
							}
						}
					}		
				}
	
			}
		} else ("x","x")

}

/* Take phased site from parents and use Homozygous SNPs in Child to drop phase down*/

def childPhase(curPhase: Tuple2[String,String], child: String): String ={
	if (child != "."){
	if ((curPhase._1 != "x") && (child == "0/0" || child == "1/1")){
		if (curPhase._1 == child(0).toString || curPhase._1 == child(1).toString) "S" else "D"
		} else {
		"U"
		}
	} else "U"

}


def getGT(data: String) : String = {
	if (data != "."){
		val tmp = data.split(":")
		if (tmp.size > GQ && tmp.size > DP && tmp(DP) != "."){ if (tmp(0) == "./." || ( tmp(GQ).toInt <= 35 && tmp(DP).toInt <= 16)) "." else tmp(0)} else "."
	} else {
		"."
	}
}

/*
*	Recombination probability function returns Tuple2(ProbDam, ProbSire)
*/

def getRecP(prevDP: Double, prevSP: Double, phase: String) : Tuple2[Double,Double] = {
val baseDamP, baseSireP = 0.49

val cD = if (phase == "D") (prevDP + (baseDamP * prevSP)) else (prevDP - (baseDamP * prevDP))
val cS = if (phase == "S") (prevSP + (baseSireP * prevDP)) else (prevSP - (baseSireP * prevSP))

(cD,cS)
}


def main (args: Array[String]) : Unit = {
import net.sf.samtools.util._

if 

val baseDamP, baseSireP = 0.49
var prevDamP, prevSireP = baseDamP
var curDamP, curSireP = 0.0


val ped = new BufferedReader(new FileReader("/scratch/aeonsim/vcfs/ref-peds/CRV-pedigree-800K.ped"))
//val GTs = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("/home/projects/bos_taurus/damona/vcfs/GATK-HC-10-Dec-2014-Mosaic-hunting.sorted.dbsnp.renamed.vcf.gz"))))
val GTs = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(new FileInputStream("/home/aeonsim/quickPhase/Damona-Mendelian-Consistant-Jul-24.vcf.gz"))))

val familyStruc = new BufferedWriter(new FileWriter(args(0) + ".txt"))
val targets = Array("1:20685483","1:20686441","1:45112139","1:51910329","1:72951069","1:97016955","1:126381931","2:75906827","2:79427071","2:112292340","2:112504415","2:134281292","3:43959715","3:53422946","3:75872618","4:113354881","5:36704515","5:41979803","5:41979804","5:120892467","6:21950345","6:35745447","6:112271783","7:21700920","7:49632157","8:10605678","8:31372887","8:40252470","9:43952072","9:86942337","9:96048216","10:491674","10:102808584","11:35306013","11:35306025","11:82969272","11:83075673","11:89422735","12:87938762","13:21363737","13:26725465","14:21325160","14:41813012","14:74925688","14:83709794","16:44658436","16:50384297","16:51332968","17:10447474","17:15961079","17:41784655","17:59076286","17:70143316","18:6505926","19:21113066","19:29060137","20:344573","20:64353624","21:64806901","21:70782567","23:13861856","23:28578136","23:29507176","23:30528411","23:52456938","25:40438599","26:10183274","26:36774923","29:3416251","X:23460318")

//val trios = Array("NL288458773","NL368790870","NL267363937","NL286575328","NL396647605","NL727555328")
//var families: List[Tuple4[String, String, String, List[String]]] = Nil

//val trios_out = trios.map(s => new BufferedWriter(new FileWriter(s + ".tab")))

/* (Sire, Dam, Child, List[Kids]) 
*   ped Columns 5 (ID), 7 (Sire) & 8 (Dam)
*/

val ID2COLUMN = new HashMap[String, Int]
var prev = new HashMap[String,String]
/*val fam1 = ("NL137409985","NL141317243","NL288458773",List("NL344991226","NL351142406","NL420753074","NL422595517","NL422743921","NL422977809","NL422992547","NL422992631","NL423125012","NL424291840","NL425596544","NL429353455","NL429486883","NL431320788","NL443684324","NL466227919","NL467115305","NL477342393","NL480773638"))*/
var currentLine = GTs.readLine.split("\t")

while (currentLine(0).apply(1) == '#') currentLine = GTs.readLine.split("\t")

for (i <- 0 to (currentLine.size -1 )){
ID2COLUMN += currentLine(i) -> i
}

val genos = ID2COLUMN.keys.toArray

var pop: List[Array[String]] = Nil
var pedFile = new HashMap[String, Array[String]]

while (ped.ready){
val temp = ped.readLine.split("\t") 
pop = temp :: pop
pedFile += temp(1) -> temp
}

if (args.size == 0) System.exit(1)
val cfam = args(0)

val anscestors = itPed(pedFile, cfam, 3)

//for (cfam <- trios){
val fam1 = (pedFile(cfam)(2),pedFile(cfam)(3),cfam,anscestors,findChildren(pop,genos,cfam))
//}
familyStruc.write("Sire\tSGS\tSGD\tDam\tDGS\tDGD\tPRoband\tAncestors\tChildren 1-> 30\n")
familyStruc.write(s"${fam1._1}")
if (pedFile.contains(fam1._1)) familyStruc.write(s"\t${pedFile(fam1._1)(2)}\t${pedFile(fam1._1)(3)}") else familyStruc.write(s"\t.\t.")
familyStruc.write(s"\t${fam1._2}")
if (pedFile.contains(fam1._2)) familyStruc.write(s"\t${pedFile(fam1._2)(2)}\t${pedFile(fam1._2)(3)}") else familyStruc.write(s"\t.\t.")
//familyStruc.write(s"${fam1._1}\t${fam1._2}\t${fam1._3}")
familyStruc.write(s"\t${fam1._3}")
//fam1._4.foreach(s => familyStruc.write("\t" + s))
fam1._5.foreach(s => familyStruc.write("\t" + s))
familyStruc.write("\n")


for(kid <- fam1._5){
prev += kid -> "."
}

var prevPos = 0
var outline = false
val outPos = new BufferedWriter(new FileWriter(cfam + ".phased.targets.txt"))
var targPos = 0

var rawPhase = new HashMap[String,List[Tuple3[String,Int,String]]]

print("CHROM\tSTART\tEND\tID" )
outPos.write("Target")
for(kid <- fam1._5){
outPos.write("\t" + kid)
print(s"\t${kid}\tD-${kid}\tS-${kid}")
rawPhase += kid -> Nil
}

outPos.write("\n")
print("\n")
var prevChrm = "0"

var prevSirePhases = new HashMap[String,Double]
var prevDamPhases = new HashMap[String,Double]
var curPhases = new HashMap[String,Double]
var prevPhase = new HashMap[String,Tuple2[String,String]]

for (kid <- fam1._5){
prevDamPhases += kid -> baseDamP
prevSirePhases += kid -> baseSireP
prevPhase += kid -> ("","")
//curDamPhases += kid -> 0.0
//curSirePhases += kid -> 0.0
}


while (GTs.ready){
currentLine = GTs.readLine.split("\t")

if (GQ == 0) GQ = currentLine(8).split(":").indexOf("GQ")
if (DP == 0) DP = currentLine(8).split(":").indexOf("DP")
val MQval = currentLine(7).split(";").indexOf("MQ")
var MQ = if (MQval > 0) currentLine(7).split(";")(MQval).toFloat else 99.0

if (currentLine(0) != prevChrm){
for (kid <- prev.keys){
prev(kid) = "."
prevDamPhases(kid) = baseDamP
prevSirePhases(kid) = baseSireP
curSireP = baseSireP
curDamP = baseDamP
}
prevChrm = currentLine(0)
prevPos = 0
}

val curPhase = if ( MQ >= 30 && currentLine(3).size == 1 && currentLine(4).size == 1) phase(getGT(currentLine(ID2COLUMN(fam1._3))),getGT(currentLine(ID2COLUMN(fam1._1))),getGT(currentLine(ID2COLUMN(fam1._2)))) else ("x","x")

/*
for (denovo <- targets){
val cand = denovo.split(":")
if (cand(0) == currentLine(1) && cand(1).toInt < currentLine(2).toInt && cand(1).toInt > prevPos){
if (outline) outPos.write(cand(0) + ":" + cand(1) + "\n") else outPos.write(cand(0) + ":" + cand(1))
outline = true
targPos = cand(1).toInt
}
} //targets
*/

if (curPhase._1 != "x"){

print(currentLine(0) + "\t" + currentLine(1) + "\t" + currentLine(1) + "\t" + currentLine(2))

for (kid <- fam1._5){

prevDamP = prevDamPhases(kid)
prevSireP = prevSirePhases(kid)


val kidphase : String = if (ID2COLUMN.contains(pedFile(kid)(2)) && ID2COLUMN.contains(pedFile(kid)(3))) {
val chPhase = phase(getGT(currentLine(ID2COLUMN(kid))),getGT(currentLine(ID2COLUMN(pedFile(kid)(2)))),getGT(currentLine(ID2COLUMN(pedFile(kid)(3)))))
if (chPhase == ("x","x")){
childPhase(curPhase,getGT(currentLine(ID2COLUMN(kid))))
} else{
familyStruc.write(currentLine(1) + "\t" + chPhase + "\n")
//Work out Which parent and what phase
if(fam1._3 == pedFile(kid)(2)){
//Sire
if (chPhase._1 == curPhase._1) "S" else if (chPhase._1 == curPhase._2) "D" else "U"
} else {
if (fam1._3 == pedFile(kid)(2)){
//Dam
if (chPhase._2 == curPhase._1) "S" else if (chPhase._2 == curPhase._2) "D" else "U"
} else "U"
}

}
} else {
childPhase(curPhase,getGT(currentLine(ID2COLUMN(kid))))
}



if (kidphase == "U") {
	print(s"\tU\t${prevDamP}\t${prevSireP}") 

	rawPhase(kid) = (currentLine(0),currentLine(1).toInt,"U") :: rawPhase(kid) 
} else {
	//def getRecP(prevDP: Double, prevSP: Double, phase: String) : Tuple2[Double,Double]
	val probs = getRecP(prevDamP,prevSireP,kidphase)
	curDamP = probs._1
	curSireP = probs._2
	
//if ((cD > cS && cD < prevSP) || (cS > cD && cS < prevDP)) true else false
	if (kidphase != prev(kid) && ((curDamP > curSireP && curDamP < prevSirePhases(kid)) || (curSireP > curDamP && curSireP < prevDamPhases(kid)))) {
		if (prevPhase(kid)._2 == kidphase ) {
		System.err.println(currentLine(0) + "\t" + prevPhase(kid)._1 + "\tBREAKPOINT\t" + kid + "\t:\t" + prev(kid) + "\t=>\t" + kidphase + s"\tDam:${prevDamPhases(kid)} -> ${curDamP}\t${prevSirePhases(kid)} -> ${curSireP}\t")
		} else{
		System.err.println(currentLine(0) + "\t" + currentLine(1) + "\tBREAKPOINT\t" + kid + "\t:\t" + prev(kid) + "\t=>\t" + kidphase + s"\tDam:${prevDamPhases(kid)} -> ${curDamP}\t${prevSirePhases(kid)} -> ${curSireP}\t")
		}
	}
	prevPhase(kid) = (currentLine(1),kidphase)
	rawPhase(kid) = (currentLine(0),currentLine(1).toInt,kidphase) :: rawPhase(kid) 
	print("\t" + kidphase + "\t" + curDamP + "\t" + curSireP)
	prev(kid) = kidphase
	prevSirePhases(kid) = curSireP
	prevDamPhases(kid) = curDamP
}

} //for kids in family

print("\n")

} 

}//Not x

prevPos = currentLine(1).toInt
outline = false
//while ready

//Input Chrom, Pos, Phase


var phaseTracking = new HashMap[String, phaseTracker]
var phaseBlock = new HashMap[String,List[Tuple5[String,Int,Int,String,Int]]]
//val chroms = Array("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","X")
val chroms = Array("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26","chr27","chr28","chr29","chrX")

for(child <- fam1._5){
	System.err.println(child + " " + rawPhase(child).size)
			var blocks = new HashMap[String,Array[Tuple5[String, Int, Int, String,Int]]]
			//var phased : List[Tuple3[String,Int,String]] = Nil 
			
			if(rawPhase(child).size > 0) {
						
			val childOrigin = rawPhase(child).reverse.toArray
			System.err.println(rawPhase(child).reverse(0) + " " + rawPhase(child).reverse(1))
			//val childOrigin = rawPhase(child).toArray
			
			var count, score = 0			
			var chrom = childOrigin(count)._1
			var parent = childOrigin(count)._3
			var start = childOrigin(count)._2
			var end = childOrigin(count)._2
			
			if (!phaseBlock.contains(child) ) phaseBlock += child -> Nil
			
			while (count < (childOrigin.size -1)){
				count +=1
				if (parent != childOrigin(count)._3 || chrom != childOrigin(count)._1){
					phaseBlock(child) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child)
					chrom = childOrigin(count)._1
					parent = childOrigin(count)._3
					start = childOrigin(count)._2
					end = childOrigin(count)._2
					score = 1
				} else {
					end = childOrigin(count)._2
					score += 1
				}				
			}
			phaseBlock(child) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child)
			
			for (ch <- chroms){
				blocks += ch -> phaseBlock(child).filter(blk => blk._1 == ch).reverse.toArray
			}		
			phaseTracking += child -> new phaseTracker(blocks)
		} else {
			for (ch <- chroms){
				blocks += ch -> Array()
			}
			phaseTracking += child -> new phaseTracker(blocks)
		}
		}

for (denovo <-  targets){
outPos.write(denovo)
val cand = denovo.split(":")
for (kid <- fam1._5){
outPos.write("\t" + (if (phaseTracking(kid).getPhase(cand(0),cand(1).toInt) == "BOTH") phaseTracking(kid).getNearestBlock(cand(0),cand(1).toInt) else phaseTracking(kid).getPhase(cand(0),cand(1).toInt)))
}
outPos.write("\n")
}



familyStruc.close
outPos.close
}//eMain

}