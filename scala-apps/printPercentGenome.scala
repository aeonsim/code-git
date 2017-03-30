import org.apache.commons.io.FileUtils._
import java.io._

val dataSet = org.apache.commons.io.FileUtils.listFiles(new File("."),Array("depths.txt"),true).iterator

		println("Proband\tAutosome%\tchrX%")

		while (dataSet.hasNext){

			val F = dataSet.next
			val fID = F.toString.split("/").last.split("_").apply(1).split("-").apply(0)
			val tmp = new BufferedReader(new FileReader(F))

			var aG, aB, xG, xB = 0.0

			aG = tmp.readLine.split("\t")(1).toDouble
			aB = tmp.readLine.split("\t")(1).toDouble
			xG = tmp.readLine.split("\t")(1).toDouble
			xB = tmp.readLine.split("\t")(1).toDouble
			tmp.close

			println(s"${fID}\t${aG}\t${xG}")

		}