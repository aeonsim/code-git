import java.io._

//Targeted concatenated model files, if de novos are named call as filter <file> <0.xx> named

val in = new BufferedReader(new FileReader (new File(args(0))))

val filterLevel = args(1).toDouble

while(in.ready){
	val line = in.readLine().split("\t")
	print(s"${line(0)}\t")
	if (args(2) == "named"){
		var count = 2
		while(count <= line.size -1){
			if (line(count).toDouble < filterLevel) print(s"${line(count -1)}\t${line(count)}\t")
			count += 2
		}
		print("\n")
	} else {
		var count = 1
		while(count <= line.size -1){
			if (line(count).toDouble < filterLevel) print(s"${line(count)}\t")
			count += 1
		}
		print("\n")
	}
}

in.close