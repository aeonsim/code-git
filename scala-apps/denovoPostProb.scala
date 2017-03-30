 

import java.io._
import scala.collection.mutable.HashMap

val inVcf = new BufferedReader(new FileReader(args(0)))
//val keepIn = new BufferedReader(new FileReader("toKeep.txt"))
val pedIn = new BufferedReader(new FileReader("fb.ped"))
val outSub = new BufferedWriter(new FileWriter("freebayes_DPP.tab"))


def rProb(prob: Double) : Double = {
  scala.math.pow(10,prob)
}


  def denovoPostProb(proPL: Array[Double], sirePL: Array[Double], damPL: Array[Double]): Tuple2[Double, Double] = {
    val mendelian3G = Array(Array(2, 2, 2), Array(2, 2, 1), Array(2, 1, 2), Array(2, 1, 1), Array(1, 2, 1), Array(1, 2, 0), Array(1, 1, 2), Array(1, 1, 1), Array(1, 1, 0), Array(1, 0, 2), Array(1, 0, 1), Array(0, 1, 1), Array(0, 1, 0), Array(0, 0, 1), Array(0, 0, 0))
    val nonmend3G = Array(Array(0, 0, 2), Array(0, 2, 0), Array(0, 2, 2), Array(0, 2, 1), Array(0, 1, 2), Array(2, 0, 0), Array(2, 0, 2), Array(2, 0, 1), Array(2, 2, 0), Array(2, 1, 0), Array(1, 0, 0), Array(1, 2, 2))
    val L_R = (rProb(proPL(0)) * rProb(sirePL(0)) * rProb(damPL(0)))
    val P_denovo = 4E-8
    val P_mendelian = 1E-3
    val P_R = (1.0 - P_denovo - P_mendelian)
    var L_mendelian, L_denovo = 0.0

    for (mend <- mendelian3G) {
      L_mendelian += (rProb(proPL(mend(0))) * rProb(sirePL(mend(1))) * rProb(damPL(mend(2))))
    }
    for (not <- nonmend3G) {
      L_denovo += (rProb(proPL(not(0))) * rProb(sirePL(not(1))) * rProb(damPL(not(2))))
    }

    val mendelian_posterior = ((L_mendelian * P_mendelian) / ((L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_R * P_R)))
    val denovo_posterior = ((L_denovo * P_denovo) / ((L_mendelian * P_mendelian) + (L_denovo * P_denovo) + (L_R * P_R)))

    (mendelian_posterior, denovo_posterior)
  }




var currentLine = inVcf.readLine.split("\t")

var keep : List[String] = Nil
var keepPos : List[Int] = Nil

/*
while(keepIn.ready){
        val tmp = keepIn.readLine
        keep = tmp :: keep
}
keepIn.close
*/

while (currentLine(0)(1) == '#') currentLine = inVcf.readLine.split("\t")

outSub.write(s"Index\tChrom\tPos\tRef\tAlt")

var idPos = new HashMap[String,Int]

for (i <- 9 to (currentLine.size - 1)){
        //outSub.write(s"\t${currentLine(i)}")
        idPos += currentLine(i) -> i
        keepPos = i :: keepPos
}


var ped : List[Array[String]] = Nil

while (pedIn.ready){
  val tmp = pedIn.readLine.split("\t")
  ped = tmp :: ped
}

for (p <- ped){
  outSub.write(s"\t${p(0)}")
}
outSub.write("\n")

while (inVcf.ready){
  currentLine = inVcf.readLine.split("\t")
  val posGL = currentLine(8).split(":").indexOf("GL")
  outSub.write(s"${currentLine(0)}:${currentLine(1)}\t${currentLine(0)}\t${currentLine(1)}\t${currentLine(3)}\t${currentLine(4)}")
  for (p <- ped){
  	if (currentLine(idPos(p(0))).split(":").size >= posGL && currentLine(idPos(p(1))).split(":").size >= posGL && currentLine(idPos(p(2))).split(":").size >= posGL){
  		if (currentLine(idPos(p(0))).split(":")(posGL) != "." && currentLine(idPos(p(1))).split(":")(posGL) != "." && currentLine(idPos(p(2))).split(":")(posGL) != "."){
		    val pro = currentLine(idPos(p(0))).split(":")(posGL).split(",").map(_.toDouble)
		    val sire = currentLine(idPos(p(1))).split(":")(posGL).split(",").map(_.toDouble)
		    val dam = currentLine(idPos(p(2))).split(":")(posGL).split(",").map(_.toDouble)

		    val dpp = if (pro.size == 3) denovoPostProb(pro,sire,dam) else (0.0,0.0)
		    outSub.write(s"\t${dpp._2.toFloat}")
		    } else {
		    	outSub.write("\t0.0")
		    }
	} else outSub.write("\t0.0")
  }
  outSub.write("\n")
}

outSub.close