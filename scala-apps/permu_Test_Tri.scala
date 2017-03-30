import java.io._

val reps = 1000000
val tris = new BufferedReader(new FileReader(new File("tris.txt")))

var triNuc : List[String] = Nil
while (tris.ready) triNuc = tris.readLine :: triNuc
tris.close
var allTris = triNuc.toSet
var results : List[Int] = Nil

for (k <- 0 to reps ){
	val data = scala.util.Random.shuffle(List.fill[String](145)("N") ++ List.fill[String](86)("M"))
	var TCT = 0
	val iter = triNuc.toIterator

	var mos, nm : List[String] = Nil
	for (i <- data){
		if (i == "N") nm = iter.next :: nm else mos = iter.next :: mos
	}

	for (i <- mos){
		if (i == "TCT>TAT" || i == "AGA>ATA") TCT += 1
	}

	results = TCT :: results

}


var sig = 0

for (j <- results) if (j == 9) sig += 1

println("Prob of 9 TCT mos is: " + sig / reps.toDouble)