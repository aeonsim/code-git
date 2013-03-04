import java.io._

object calls2bed {

def main(args : Array[String]) : Unit = {
val input = new BufferedReader(new FileReader(args(0)))
val output = new BufferedWriter(new FileWriter(args(0) + ".bed"))

while (input.ready){
val cline = input.readLine.split("\t")
val pos = cline(1).split(":")
output.write(pos(0) + "\t" + pos(1).split("-")(0) + "\t" + pos(1).split("-")(1) + "\t" + cline(0) + "\t" + cline(4) + "\n")
}
input.close
output.close
}

}