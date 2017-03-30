import java.io._
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

val in = new BufferedReader(new FileReader(new File(args(0))))
val out = new BufferedWriter(new FileWriter(new File(args(0).split("\\.").apply(0) + "_ARform.tab")))

var line = in.readLine.split("\t")
var headers = new HashMap[Int,String]
var headersRev = new HashMap[String,Int]

var phaseTrack = new HashMap[String,String]

var haveStart = false
var dataStart = 0
var order : List[String] = Nil
var labelSet = new HashSet[String]

var seqPCounts, seqMCounts, valPCounts , valMCounts = new HashMap[String,String]

for (i <- 0 to (line.size - 1)){
	headers += i -> line(i)
	headersRev += line(i) -> i
	if (!haveStart && line(i).contains("PR")) {
		dataStart = i
		haveStart = true
	}
}

def calcAR(AR: String) : String = {
	if (AR != '-' && AR.contains("|")){
			( AR.substring(AR.indexOf("|") + 1, AR.size).toFloat / (AR.substring(AR.indexOf("|") + 1, AR.size).toFloat + AR.substring(0,AR.indexOf("|")).toFloat)).toString
		} else {
			"-"
		}
}

def aForm(dat: String, head: String, id: String, name: String): Unit = {
	out.write(s"${head}${id}\t${name}\t${calcAR(dat)}\t-\n")
}

def aRel(dat: String, head: String, id: String, name: String,phase : String): Unit = {
	var lar = calcAR(dat)
	if (lar != "-" && lar.toFloat <= 0.01) phase match {
		case x if x == "M" => lar = "-0.01"
		case x if x == "P" => lar = "0.01"
		case _ => 
	}
	if (phase != "-") out.write(s"${head}${id + " (" + phase + ")"}\t${name}\t${lar}\t${phase}\n")
}

out.write("Class\tMutation\tPosition\tType\tGrouping\tName\tAR\tPhase\n")
var other :List[Int] = Nil

while (in.ready){
	line = in.readLine.split("\t")
	val descriptor = s"${line(4)}\t${line(0).substring(3,line(0).size)}\t${line(2)}:${line(3)}\t${line(1).toUpperCase.head}\t"

	labelSet += line(0).substring(3,line(0).size)

	order = line(0).substring(3,line(0).size) :: order
	val mutClass = line(4)
	for (i <- dataStart to (line.size -1)){
		headers(i) match{
			case x if x.toUpperCase.contains("SIR") => aForm(line(i),descriptor,"Sire",x) //if (List(2,3).contains(mutClass)) aForm(line(i),descriptor,"Sire",x)
			case x if x.toUpperCase.contains("DAM") => aForm(line(i),descriptor,"Dam",x) //if (List(4,5).contains(mutClass)) aForm(line(i),descriptor,"Dam",x)
			case x if x.toUpperCase.contains("PGS") => if (mutClass == "3") aForm(line(i),descriptor,"PGS",x)
			case x if x.toUpperCase.contains("PGD") => if (mutClass == "3") aForm(line(i),descriptor,"PGD",x)
			case x if x.toUpperCase.contains("MGS") => if (mutClass == "5") aForm(line(i),descriptor,"MGS",x)
			case x if x.toUpperCase.contains("MGD") => if (mutClass == "5") aForm(line(i),descriptor,"MGD",x)
			case x if x.toUpperCase.contains("PRO") => aForm(line(i),descriptor,"Proband",x)
			case _ => other = i :: other
		}
	}

	if (List("1","2","4").contains(mutClass)){
			val goPM = line(headersRev("GO Paternal Mut"))
			val goPH = line(headersRev("GO Paternal Hap"))
			val goMM = line(headersRev("GO Maternal Mut"))
			val goMH = line(headersRev("GO Maternal Hap"))

			if (line(1).toUpperCase.head == 'S') {if (goPM.toInt >= 1) phaseTrack += line(0).substring(3,line(0).size) -> "P" else phaseTrack += line(0).substring(3,line(0).size) -> "M"}

			if (goPM != "-"){
				if (line(1).toUpperCase.head == 'S') seqPCounts += line(0).substring(3,line(0).size) -> (goPH.toInt - goPM.toInt).toString
				if (line(1).toUpperCase.head == 'S') seqMCounts += line(0).substring(3,line(0).size) -> (goMH.toInt - goMM.toInt).toString

				if (line(1).toUpperCase.head == 'V') valPCounts += line(0).substring(3,line(0).size) -> (goPH.toInt - goPM.toInt).toString
				if (line(1).toUpperCase.head == 'V') valMCounts += line(0).substring(3,line(0).size) -> (goMH.toInt - goMM.toInt).toString
			} else {
				if (line(1).toUpperCase.head == 'S') {
					seqPCounts += line(0).substring(3,line(0).size) -> "-"
					seqMCounts += line(0).substring(3,line(0).size) -> "-"
				}
				if (line(1).toUpperCase.head == 'V') {
					valPCounts += line(0).substring(3,line(0).size) -> "-"
					valMCounts += line(0).substring(3,line(0).size) -> "-"
				}

			}


		} else {

			val goPM = line(headersRev("HS Paternal Mut"))
			val goPH = line(headersRev("HS Paternal Hap"))
			val goMM = line(headersRev("HS Maternal Mut"))
			val goMH = line(headersRev("HS Maternal Hap"))

			if (line(1).toUpperCase.head == 'S') {if (goPM.toInt >= 1) phaseTrack += line(0).substring(3,line(0).size) -> "P" else phaseTrack += line(0).substring(3,line(0).size) -> "M"}

			if (goPM != "-"){
				if (line(1).toUpperCase.head == 'S') seqPCounts += line(0).substring(3,line(0).size) -> (goPH.toInt - goPM.toInt).toString
				if (line(1).toUpperCase.head == 'S') seqMCounts += line(0).substring(3,line(0).size) -> (goMH.toInt - goMM.toInt).toString

				if (line(1).toUpperCase.head == 'V') valPCounts += line(0).substring(3,line(0).size) -> (goPH.toInt - goPM.toInt).toString
				if (line(1).toUpperCase.head == 'V') valMCounts += line(0).substring(3,line(0).size) -> (goMH.toInt - goMM.toInt).toString
			} else {
				if (line(1).toUpperCase.head == 'S') {
					seqPCounts += line(0).substring(3,line(0).size) -> "-"
					seqMCounts += line(0).substring(3,line(0).size) -> "-"
				}
				if (line(1).toUpperCase.head == 'V') {
					valPCounts += line(0).substring(3,line(0).size) -> "-"
					valMCounts += line(0).substring(3,line(0).size) -> "-"
				}

			}

		}



	var count = other.sorted.head
	val end = other.sorted.last

	while (count < end ){
		headers(count) match {
			case x if x.toUpperCase.contains("GO") => aRel(line(count),descriptor,"Offspring",x,line(count + 1))
			case x if x.toUpperCase.contains("PHS") => if (mutClass == "3") aRel(line(count),descriptor,"Pat Halfsib",x,line(count + 1))
			case x if x.toUpperCase.contains("MHS") => if (mutClass == "5") aRel(line(count),descriptor,"Mat Halfsib",x,line(count + 1))
			case _ => println(headers(count))
		}

		count += 2
	}

}

order.reverse.foreach(s => print( "\"" + s + "\","))
print("\n\n")

for (i <- labelSet) {
	if (phaseTrack(i) == "P") {
			print("`" + i + "` = \"" + i + "\\n" + seqPCounts(i) + "   " + valPCounts(i) + "\\n" + seqMCounts(i) + "   " + valMCounts(i) + "\",") 
		} else { 
			print("`" + i + "` = \"" + i + "\\n" + seqMCounts(i) + "   " + valMCounts(i) + "\\n" + seqPCounts(i) + "   " + valPCounts(i) + "\",")
		}
}
print("\n\n")

out.close