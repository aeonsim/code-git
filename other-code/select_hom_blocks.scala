import java.io._

var start, ff, stop, last, cnt, het, hom = 0
var prev = 0.0
var cChr = ""

val inputs = Array("animal_10841297.bedgraph",
"animal_11056944.bedgraph",
"animal_11955102.bedgraph",
"animal_13192605.bedgraph",
"animal_14415740.bedgraph",
"animal_14915827.bedgraph",
"animal_15525685.bedgraph",
"animal_15656114.bedgraph",
"animal_16109405.bedgraph",
"animal_16120004.bedgraph",
"animal_16703698.bedgraph",
"animal_17034899.bedgraph",
"animal_17074197.bedgraph",
"animal_17108826.bedgraph",
"animal_17144783.bedgraph",
"animal_17144784.bedgraph",
"animal_17145716.bedgraph",
"animal_17396939.bedgraph",
"animal_18033109.bedgraph",
"animal_18093560.bedgraph",
"animal_18143602.bedgraph",
"animal_18186561.bedgraph",
"animal_20976747.bedgraph",
"animal_22829340.bedgraph",
"animal_6990172.bedgraph")

for (ans <- inputs){
val auto_prob = new BufferedReader(new FileReader(ans))
val out = new BufferedWriter(new FileWriter(ans.split("\\.").apply(0) + "homozygous_blocks.txt"))

def endblock (chr: String, end: Int) : Unit = {
out.write(chr + "\t" + start + "\t" + end + "\t" + het + "\t" + hom + "\t" + (end - start) + "\n")
start = 0
het = 0
hom = 0
prev = 0.0
cChr = ""
cnt = 0
last = 0
}
val cur = auto_prob.readLine.split("\t")

while (auto_prob.ready){
val cur = auto_prob.readLine.split("\t")

if ((cur(2).toInt - start) < 0){
endblock(cChr,last)
}

if (start != 0){
if (cur(3).toFloat >= 0.99){
last = cur(2).toInt
prev = cur(3).toFloat
hom += 1
} else {
if (het == 0){
prev = cur(3).toFloat
het += 1
ff = cur(2).toInt
} else {
het += 1
}
if (het >= 3){
endblock(cur(0),cur(2).toInt)
}
}
}else{
if (cur(3).toFloat >= 0.999) {
println("Starting new Region \n")
start = cur(2).toInt
prev = cur(3).toFloat
hom += 1
cChr = cur(0)
}
}
}
auto_prob.close
out.close

}
