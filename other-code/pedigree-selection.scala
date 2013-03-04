import scala.io._
import java.io._
import scala.collection.mutable.HashMap

object trios {

def main(args:Array[String]): Unit = {

if (args.size == 0){
println("trio <min children> <filepath>")
System.exit(0)
}

val rawData = new BufferedReader(new FileReader(args(1)))

val trios = new HashMap[String, Int]
//Animal -> Sire Dam inbreed NA-HF Friesian16ths
val data = new HashMap[String,Tuple7[String,String,String,String,String,String,String]]
val children = new HashMap[String,Int]

/*
LIC_ANIMAL_KEY 0,name 1,sire_lic_animal_key 2,dam_lic_animal_key 3,sex_code 4,birth_id_year 5,
friesian_breed_16ths 6,jersey_breed_16ths 7,inbreed_coefficient 8,north_american_hf_pct 9
*/
println(rawData.readLine)

while (rawData.ready) {
val line = rawData.readLine.split(",")
//Animal -> Sire Dam inbreed NA-HF HF16ths sex dob
if (line(6).toInt >= 14){
data += line(0) -> (line(2),line(3),line(8),line(9),line(6),line(4),line(5))
if (children.contains(line(2))){
children(line(2)) += 1
} else {
children += line(2) -> 1
}
if (children.contains(line(3))){
children(line(3)) += 1
} else {
children += line(3) -> 1
}
}
}

for (f1 <- data.keys){
if(children.contains(f1)){
//if ((children(f1) >= 5) & (data.contains(data(f1)._1)) & (data.contains(data(f1)._2))){
if ((children(f1) >= args(0).toInt) & (data.contains(data(f1)._1)) & (data.contains(data(f1)._2))){
trios += f1 -> 1
}
}
}

for (f1 <- trios.keys){
if(trios.contains(data(f1)._1)){
trios(f1) += 1
trios(data(f1)._1) += 1
}
if(trios.contains(data(f1)._2)){
trios(f1) += 1
trios(data(f1)._2) += 1
}
}

val outTrios = new BufferedWriter(new FileWriter("NZ_HF_Trios.txt"))

outTrios.write("F1\tChildren\tSire\tDAM\tSex\tDOB\tHF16ths\tInbCo\tNAHF%\tPedigreeSize\n")

for (f1 <- trios.keys){
val dtls = data(f1)
outTrios.write(f1 + "\t" + children(f1) + "\t" + dtls._1  + "\t" + dtls._2  + "\t" + dtls._6  + "\t" + dtls._7  + "\t" + dtls._5  + "\t" + dtls._3  + "\t" + dtls._4 + "\t" + trios(f1) + "\n")
}
outTrios.close

//Check Multi Generations

val outDetail = new BufferedWriter(new FileWriter("NZ_HF_Trios-details.txt"))

for (f1 <- trios.keys){
def check(anml:String) : Int = {
var cgen = 0
if (trios.contains(anml)){
cgen += check(data(anml)._1)
cgen += check(data(anml)._2)
}
cgen + 1
}
var gens = check(f1)
outDetail.write(f1 + "\t" + gens + "\n")
}


for (f1 <- trios.keys){
def check(anml:String) : String = {
var cgen = ""
if (trios.contains(anml)){
cgen = cgen + "," + check(data(anml)._1)
cgen = cgen + "," + check(data(anml)._2)
}
cgen + "," +anml
}
var gens = check(f1)
outDetail.write(f1 + "\t" + gens)
}

outDetail.close
val outCyto = new BufferedWriter(new FileWriter("NZ_HF_Trios-Cyto.txt"))

for (f1 <- trios.keys){

def check2(F0:String,pro:String) : Unit = {
if (trios.contains(F0)){
if (pro != "") outCyto.write(pro + "\t" + F0 + "\n")
check2(data(F0)._1,F0)
check2(data(F0)._2,F0)
}else{
if (pro != "") outCyto.write(pro + "\t" + F0 + "\n")
}
}
check2(f1,"")
}
outCyto.close
}
}
