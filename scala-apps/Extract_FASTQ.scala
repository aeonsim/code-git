import java.io._
import scala.collection.mutable.HashSet

System.err.println("scala app.scala <targetReadNames> <targetFastq>")

val inTargets = new BufferedReader(new FileReader(new File(args(0))))
val inFastQ = new BufferedReader(new FileReader(new File(args(1))))

//val inTargets = new BufferedReader(new FileReader(new File(~/readNames.targ)))
//val inFastQ = new BufferedReader(new FileReader(new File(17135_filtered_subreads.fastq)))

var targNames = new HashSet[String]

var cline = ""

while (inTargets.ready){
	cline = inTargets.readLine
	targNames += ("@" + cline)
}


cline = " "

while(inFastQ.ready){
	if (cline(0) == '@'){
		if (targNames.contains(cline)){
			println(cline)
			cline = inFastQ.readLine
			while (cline(0) != '@') {
				println(cline)
				cline = inFastQ.readLine
			}
		} else {
			cline = inFastQ.readLine
		}
	} else {
		cline = inFastQ.readLine
	}
}