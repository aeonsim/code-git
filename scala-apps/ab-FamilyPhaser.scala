
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


object abPhaser {

import java.io._
import scala.collection.mutable.HashMap

/* Recursive function to find all parents & ancestors of nominated individual
* returns list of those present in both the PED & VCF files.
*/


	def itPed (pop: HashMap[String, Array[String]] , rec: String, maxIT: Int): List[String] = {
		if (pop.contains(rec) && maxIT >= 0){
			return rec :: itPed(pop,pop(rec).apply(7),maxIT -1) ::: itPed(pop,pop(rec).apply(8),maxIT -1)
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
			if (vcf.contains(item(5)) && ( (item(7) == parent)||(item(8) == parent) )){
				result = item(5) :: result
			}//eif

		}//efor
		result
	}//edef
	
/* Phase Code, return format is (sireAllele,damAllele)*/

def phase(indv: String, sire: String, dam: String) : Tuple2[String, String] = {
if (sire != "00" && dam != "00" && indv != "00"){
	if(sire == "AA" && dam == "BB" && (indv == "AB" || indv == "BA")){
		("A","B")
	} else {
	 	if (sire == "BB" && dam == "AA" && (indv == "AB" || indv == "BA")){
	 		("B","A")
	 		} else {
				if(sire == "AB" && (indv == "AB" || indv == "BA") && dam == "AA") {
					("B","A")
				} else {
					if(sire == "AB" && (indv == "AB" || indv == "BA") && dam == "BB") {
						("A","B")
					} else {
						if (dam == "AB" && (indv == "AB" || indv == "BA") && sire == "AA" ) {
							("A","B")
						} else {
							if (dam == "AB" && (indv == "AB" || indv == "BA") && sire == "BB") {
								("B","A")
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
	if (child != "00"){
	if ((curPhase._1 != "x") && (child == "AA" || child == "BB")){
		if (curPhase._1 == child(0).toString || curPhase._1 == child(1).toString) "S" else "D"
		} else {
		"U"
		}
	} else "U"

}


def main (args: Array[String]) : Unit = {
val ped = new BufferedReader(new FileReader("/Users/chhar0/Dropbox/PHD-WORK/VALIDATION/Validation/ValPed_20150428.txt"))
val GTs = new BufferedReader(new FileReader("/Users/chhar0/Dropbox/PHD-WORK/VALIDATION/Validation/ValGenos_20150428_sorted.txt"))
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
//val fam1 = ("NL137409985","NL141317243","NL288458773",List("NL344991226","NL351142406","NL420753074","NL422595517","NL422743921","NL422977809","NL422992547","NL422992631","NL423125012","NL424291840","NL425596544","NL429353455","NL429486883","NL431320788","NL443684324","NL466227919","NL467115305","NL477342393","NL480773638"))
var currentLine = GTs.readLine.split("\t")

for (i <- 0 to (currentLine.size -1 )){
ID2COLUMN += currentLine(i) -> i
}

val genos = ID2COLUMN.keys.toArray

var pop: List[Array[String]] = Nil
var pedFile = new HashMap[String, Array[String]]

while (ped.ready){
val temp = ped.readLine.split("\t") 
pop = temp :: pop
pedFile += temp(5) -> temp
}

if (args.size == 0) System.exit(1)
val cfam = args(0)

val anscestors = itPed(pedFile, cfam, 3)

//for (cfam <- trios){
val fam1 = (pedFile(cfam)(7),pedFile(cfam)(8),cfam,anscestors,findChildren(pop,genos,cfam))
//}
familyStruc.write("Sire\tSGS\tSGD\tDam\tDGS\tDGD\tPRoband\tAncestors\tChildren 1-> 30\n")
familyStruc.write(s"${fam1._1}")
if (pedFile.contains(fam1._1)) familyStruc.write(s"\t${pedFile(fam1._1)(7)}\t${pedFile(fam1._1)(8)}") else familyStruc.write(s"\t.\t.")
familyStruc.write(s"\t${fam1._2}")
if (pedFile.contains(fam1._2)) familyStruc.write(s"\t${pedFile(fam1._2)(7)}\t${pedFile(fam1._2)(8)}") else familyStruc.write(s"\t.\t.")
//familyStruc.write(s"${fam1._1}\t${fam1._2}\t${fam1._3}")
familyStruc.write(s"\t${fam1._3}")
//fam1._4.foreach(s => familyStruc.write("\t" + s))
fam1._5.foreach(s => familyStruc.write("\t" + s))
familyStruc.write("\n")
familyStruc.close

for(kid <- fam1._5){
prev += kid -> "."
}

var prevPos = 0
var outline = false
val outPos = new BufferedWriter(new FileWriter(cfam + ".phased.targets.txt"))
var targPos = 0

print("ID\tChr\tPos\tGT1\tGT2\t=>" )
outPos.write("Target")
for(kid <- fam1._5){
print("\t" + kid)
outPos.write("\t" + kid)
}
outPos.write("\n")
print("\n")
var prevChrm = "0"

while (GTs.ready){
currentLine = GTs.readLine.split("\t")

if (currentLine(1) != prevChrm){
for (kid <- prev.keys){
prev(kid) = "."
}
prevChrm = currentLine(1)
prevPos = 0
}

val curPhase = phase(currentLine(ID2COLUMN(fam1._3)),currentLine(ID2COLUMN(fam1._1)),currentLine(ID2COLUMN(fam1._2)))


for (denovo <- targets){
val cand = denovo.split(":")
if (cand(0) == currentLine(1) && cand(1).toInt < currentLine(2).toInt && cand(1).toInt > prevPos){
if (outline) outPos.write(cand(0) + ":" + cand(1) + "\n") else outPos.write(cand(0) + ":" + cand(1))
outline = true
targPos = cand(1).toInt
}
} //targets


if (curPhase._1 != "x"){

print(currentLine(0) + " " + currentLine(1) + " " + currentLine(2) + " " + curPhase._1 + " " + curPhase._2 + " => ")

for (kid <- fam1._5){

val kidphase = childPhase(curPhase,currentLine(ID2COLUMN(kid)))

if (kidphase != prev(kid) && kidphase != "U") {
System.err.println(currentLine(1) + "\t" + currentLine(2) + "\tBREAKPOINT\t" + kid + "\t:\t" + prev(kid) + "\t=>\t" + kidphase)

val prevDif = scala.math.abs(prevPos - targPos)
val curDif = scala.math.abs(currentLine(2).toInt - targPos)

if (outline) { 
if (prev(kid) == ".") outPos.write("\t" + kid + "\t" + kidphase) else {
if (prevDif >= curDif) outPos.write("\t" + kid + "\t" + kidphase) else outPos.write("\t" + kid + "\t" + prev(kid)) 
}
}

}

if (kidphase == "U") {

if (outline) outPos.write("\t" + prev(kid))
print("\t" + kid + "\t" + prev(kid)) 

} else {

if (outline) outPos.write("\t" + kid + "\t" + kidphase)
print("\t" + kidphase)
prev(kid) = kidphase

}

} //for kids in family

print("\n")
if (outline) outPos.write("\n")

} else {
if (outline){
for (kid <- fam1._5){
outPos.write("\t" + kid + "\t" + prev(kid))
}
outPos.write("\n")
}

}//Not x

prevPos = currentLine(2).toInt
outline = false
}//while ready

//Input Chrom, Pos, Phase

for(child <- fam1._5){
			var blocks = new HashMap[String,Array[Tuple5[String, Int, Int, String,Int]]]
			var phased : List[Tuple3[String,Int,String]] = Nil 
			
		while (in.ready){
				/* Read file into Array of tuples */
				val line = in.readLine.split("\t")
				if (line.size >= 3) phased = (line(0),line(1).toInt,line(2)) :: phased
			}
			in.close
			
			if(phased.size > 0) {
						
			val childOrigin = phased.reverse.toArray
			
			var count, score = 0			
			var chrom = childOrigin(count)._1
			var parent = childOrigin(count)._3
			var start = childOrigin(count)._2
			var end = childOrigin(count)._2
			
			if (!phaseBlock.contains(child._1) ) phaseBlock += child._1 -> Nil
			
			while (count < (childOrigin.size -1)){
				count +=1
				if (parent != childOrigin(count)._3 || chrom != childOrigin(count)._1){
					phaseBlock(child._1) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child._1)
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
			phaseBlock(child._1) = Tuple5(chrom,start,end,parent,score) :: phaseBlock(child._1)
			
			for (ch <- vcfChrs){
				blocks += ch -> phaseBlock(child._1).filter(blk => blk._5 > 2).filter(blk => blk._1 == ch).reverse.toArray
			}		
			phaseTracking += child._1 -> new phaseTracker(blocks)
		} else {
			for (ch <- vcfChrs){
				blocks += ch -> Array()
			}
			phaseTracking += child._1 -> new phaseTracker(blocks)
		}
	}

phaseBlock = new HashMap[String,List[Tuple5[String,Int,Int,String,Int]]]




outPos.close
}//eMain

}