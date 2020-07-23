suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
library(BiocFileCache)
library(curatedMetagenomicData)
library(readxl)
library(skimr)

# The combined Metadata Excel sheet should:
#   1. Have a column called "dataset_name" with a reference to the researcher or name of the study associated
#   2. Have a column named "disease" (where 0 indicates control or lack of T2D, 1 indicates prediabetes, and 2 indicates diabetes)

main <- function() {
  ## QIN ###
  # Part 1: From curatedMetagenomic Data
  qin_data_samples <- curatedMetagenomicData("QinJ_2012.metaphlan_bugs_list.stool", dryrun=FALSE)
  qin_data_samples <- pData(qin_data_samples[[1]])
  qin_data_samples <- qin_data_samples %>%
    rename(sampleID = subjectID)
  
  # Part 2, from another paper's supplementary materials
  GET("ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100064/gut_microbiome_samples_20140613.xlsx",
      write_disk(tf <- tempfile(fileext = ".xlsx")))
  qin_data_patient <- read_excel(tf, 1L)
  qin_data_patient <- qin_data_patient %>% select(Name, SRASample) %>%
    rename(sampleID = Name)
  
  temp_join <- left_join(qin_data_samples, qin_data_patient, by = "sampleID")
  
  # Part 3, from the supplementary materials of the paper
  GET("https://static-content.springer.com/esm/art%3A10.1038%2Fnature11450/MediaObjects/41586_2012_BFnature11450_MOESM372_ESM.xls",
      write_disk(tf <- tempfile(fileext = ".xls")))
  supplementary_info <- read_excel(tf, 2L)
  supplementary_info <- supplementary_info %>%
    select( "Sample ID":"LDL (mmol/L)") %>%
    rename(sampleID = "Sample ID")
  
  temp_join <- full_join(temp_join, supplementary_info, by = "sampleID")
  
  # Part 4, from the download of the raw data on the terminal
  Qin_key <- tribble(
    ~SRASample, ~SRR,
    "SRS333664", "SRR1778450",
    "SRS333665", "SRR1778451",
    "SRS333666", "SRR1778452",
    "SRS333667", "SRR1778453",
    "SRS333668", "SRR1778454",
    "SRS333669", "SRR1778455",
    "SRS333670", "SRR1778456",
    "SRS259434", "SRR341581",
    "SRS259435", "SRR341582",
    "SRS259436", "SRR341583",
    "SRS259437", "SRR341584",
    "SRS259438", "SRR341585",
    "SRS259439", "SRR341586",
    "SRS259440", "SRR341587",
    "SRS259441", "SRR341588",
    "SRS259442", "SRR341589",
    "SRS259443", "SRR341590",
    "SRS259444", "SRR341591",
    "SRS259445", "SRR341592",
    "SRS259446", "SRR341593",
    "SRS259447", "SRR341594",
    "SRS259448", "SRR341595",
    "SRS259449", "SRR341596",
    "SRS259450", "SRR341597",
    "SRS259451", "SRR341598",
    "SRS259452", "SRR341599",
    "SRS259453", "SRR341600",
    "SRS259454", "SRR341601",
    "SRS259455", "SRR341602",
    "SRS259456", "SRR341603",
    "SRS259457", "SRR341604",
    "SRS259458", "SRR341605",
    "SRS259459", "SRR341606",
    "SRS259460", "SRR341607",
    "SRS259461", "SRR341608",
    "SRS259462", "SRR341609",
    "SRS259463", "SRR341610",
    "SRS259464", "SRR341611",
    "SRS259465", "SRR341612",
    "SRS259466", "SRR341613",
    "SRS259467", "SRR341614",
    "SRS259468", "SRR341615",
    "SRS259469", "SRR341616",
    "SRS259470", "SRR341617",
    "SRS259471", "SRR341618",
    "SRS259472", "SRR341619",
    "SRS259473", "SRR341620",
    "SRS259474", "SRR341621",
    "SRS259475", "SRR341622",
    "SRS259476", "SRR341623",
    "SRS259477", "SRR341624",
    "SRS259478", "SRR341625",
    "SRS259479", "SRR341626",
    "SRS259480", "SRR341627",
    "SRS259481", "SRR341628",
    "SRS259482", "SRR341629",
    "SRS259483", "SRR341630",
    "SRS259484", "SRR341631",
    "SRS259485", "SRR341632",
    "SRS259486", "SRR341633",
    "SRS259487", "SRR341634",
    "SRS259488", "SRR341635",
    "SRS259489", "SRR341636",
    "SRS259490", "SRR341637",
    "SRS259491", "SRR341638",
    "SRS259492", "SRR341639",
    "SRS259493", "SRR341640",
    "SRS259494", "SRR341641",
    "SRS259495", "SRR341642",
    "SRS259496", "SRR341643",
    "SRS259497", "SRR341644",
    "SRS259498", "SRR341645",
    "SRS259499", "SRR341646",
    "SRS259500", "SRR341647",
    "SRS259501", "SRR341648",
    "SRS259502", "SRR341649",
    "SRS259503", "SRR341650",
    "SRS259504", "SRR341651",
    "SRS259505", "SRR341652",
    "SRS259506", "SRR341653",
    "SRS259507", "SRR341654",
    "SRS259508", "SRR341655",
    "SRS259509", "SRR341656",
    "SRS259510", "SRR341657",
    "SRS259511", "SRR341658",
    "SRS259512", "SRR341659",
    "SRS259513", "SRR341660",
    "SRS259514", "SRR341661",
    "SRS259515", "SRR341662",
    "SRS259516", "SRR341663",
    "SRS259517", "SRR341664",
    "SRS259518", "SRR341665",
    "SRS259519", "SRR341666",
    "SRS259520", "SRR341667",
    "SRS259521", "SRR341668",
    "SRS259522", "SRR341669",
    "SRS259523", "SRR341670",
    "SRS259524", "SRR341671",
    "SRS259525", "SRR341672",
    "SRS259526", "SRR341673",
    "SRS259527", "SRR341674",
    "SRS259528", "SRR341675",
    "SRS259529", "SRR341676",
    "SRS259530", "SRR341677",
    "SRS259531", "SRR341678",
    "SRS259532", "SRR341679",
    "SRS259533", "SRR341680",
    "SRS259534", "SRR341681",
    "SRS259535", "SRR341682",
    "SRS259536", "SRR341683",
    "SRS259537", "SRR341684",
    "SRS259538", "SRR341685",
    "SRS259539", "SRR341686",
    "SRS259540", "SRR341687",
    "SRS259541", "SRR341688",
    "SRS259542", "SRR341689",
    "SRS259543", "SRR341690",
    "SRS259544", "SRR341691",
    "SRS259545", "SRR341692",
    "SRS259546", "SRR341693",
    "SRS259547", "SRR341694",
    "SRS259548", "SRR341695",
    "SRS259549", "SRR341696",
    "SRS259550", "SRR341697",
    "SRS259551", "SRR341698",
    "SRS259552", "SRR341699",
    "SRS259553", "SRR341700",
    "SRS259554", "SRR341701",
    "SRS259555", "SRR341702",
    "SRS259556", "SRR341703",
    "SRS259557", "SRR341704",
    "SRS259558", "SRR341705",
    "SRS259559", "SRR341706",
    "SRS259560", "SRR341707",
    "SRS259561", "SRR341708",
    "SRS259562", "SRR341709",
    "SRS259563", "SRR341710",
    "SRS259564", "SRR341711",
    "SRS259565", "SRR341712",
    "SRS259566", "SRR341713",
    "SRS259567", "SRR341714",
    "SRS259568", "SRR341715",
    "SRS259569", "SRR341716",
    "SRS259570", "SRR341717",
    "SRS259571", "SRR341718",
    "SRS259572", "SRR341719",
    "SRS259573", "SRR341720",
    "SRS259574", "SRR341721",
    "SRS259575", "SRR341722",
    "SRS259576", "SRR341723",
    "SRS259577", "SRR341724",
    "SRS259578", "SRR341725",
    "SRS294812", "SRR413556",
    "SRS294813", "SRR413557",
    "SRS294814", "SRR413558",
    "SRS294815", "SRR413559",
    "SRS294816", "SRR413560",
    "SRS294817", "SRR413561",
    "SRS294818", "SRR413562",
    "SRS294819", "SRR413563",
    "SRS294820", "SRR413564",
    "SRS294821", "SRR413565",
    "SRS294822", "SRR413566",
    "SRS294823", "SRR413567",
    "SRS294824", "SRR413568",
    "SRS294825", "SRR413569",
    "SRS294826", "SRR413570",
    "SRS294827", "SRR413571",
    "SRS294828", "SRR413572",
    "SRS294829", "SRR413573",
    "SRS294830", "SRR413574",
    "SRS294831", "SRR413575",
    "SRS294832", "SRR413576",
    "SRS294833", "SRR413577",
    "SRS294834", "SRR413578",
    "SRS294835", "SRR413579",
    "SRS294836", "SRR413580",
    "SRS294837", "SRR413581",
    "SRS294838", "SRR413582",
    "SRS294839", "SRR413583",
    "SRS294840", "SRR413584",
    "SRS294841", "SRR413585",
    "SRS294842", "SRR413586",
    "SRS294843", "SRR413587",
    "SRS294844", "SRR413588",
    "SRS294845", "SRR413589",
    "SRS294846", "SRR413590",
    "SRS294847", "SRR413591",
    "SRS294848", "SRR413592",
    "SRS294849", "SRR413593",
    "SRS294850", "SRR413594",
    "SRS294851", "SRR413595",
    "SRS294852", "SRR413596",
    "SRS294853", "SRR413597",
    "SRS294854", "SRR413598",
    "SRS294855", "SRR413599",
    "SRS294856", "SRR413600",
    "SRS294857", "SRR413601",
    "SRS294858", "SRR413602",
    "SRS294859", "SRR413603",
    "SRS294860", "SRR413604",
    "SRS294861", "SRR413605",
    "SRS294862", "SRR413606",
    "SRS294863", "SRR413607",
    "SRS294864", "SRR413608",
    "SRS294865", "SRR413609",
    "SRS294866", "SRR413610",
    "SRS294867", "SRR413611",
    "SRS294868", "SRR413612",
    "SRS294869", "SRR413613",
    "SRS294870", "SRR413614",
    "SRS294871", "SRR413615",
    "SRS294872", "SRR413616",
    "SRS294873", "SRR413617",
    "SRS294874", "SRR413618",
    "SRS294875", "SRR413619",
    "SRS294876", "SRR413620",
    "SRS294877", "SRR413621",
    "SRS294878", "SRR413622",
    "SRS294879", "SRR413623",
    "SRS294880", "SRR413624",
    "SRS294881", "SRR413625",
    "SRS294882", "SRR413626",
    "SRS294883", "SRR413627",
    "SRS294884", "SRR413628",
    "SRS294885", "SRR413629",
    "SRS294886", "SRR413630",
    "SRS294887", "SRR413631",
    "SRS294888", "SRR413632",
    "SRS294889", "SRR413633",
    "SRS294890", "SRR413634",
    "SRS294891", "SRR413635",
    "SRS294892", "SRR413636",
    "SRS294893", "SRR413637",
    "SRS294894", "SRR413638",
    "SRS294895", "SRR413639",
    "SRS294896", "SRR413640",
    "SRS294897", "SRR413641",
    "SRS294898", "SRR413642",
    "SRS294899", "SRR413643",
    "SRS294900", "SRR413644",
    "SRS294901", "SRR413645",
    "SRS294902", "SRR413646",
    "SRS294903", "SRR413647",
    "SRS294904", "SRR413648",
    "SRS294905", "SRR413649",
    "SRS294906", "SRR413650",
    "SRS294907", "SRR413651",
    "SRS294908", "SRR413652",
    "SRS294909", "SRR413653",
    "SRS294910", "SRR413654",
    "SRS294911", "SRR413655",
    "SRS294912", "SRR413656",
    "SRS294913", "SRR413657",
    "SRS294914", "SRR413658",
    "SRS294915", "SRR413659",
    "SRS294916", "SRR413660",
    "SRS294917", "SRR413661",
    "SRS294918", "SRR413662",
    "SRS294919", "SRR413663",
    "SRS294920", "SRR413664",
    "SRS294921", "SRR413665",
    "SRS294922", "SRR413666",
    "SRS294923", "SRR413667",
    "SRS294924", "SRR413668",
    "SRS294925", "SRR413669",
    "SRS294926", "SRR413670",
    "SRS294927", "SRR413671",
    "SRS294928", "SRR413672",
    "SRS294929", "SRR413673",
    "SRS294930", "SRR413674",
    "SRS294931", "SRR413675",
    "SRS294932", "SRR413676",
    "SRS294933", "SRR413677",
    "SRS294934", "SRR413678",
    "SRS294935", "SRR413679",
    "SRS294936", "SRR413680",
    "SRS294937", "SRR413681",
    "SRS294938", "SRR413682",
    "SRS294939", "SRR413683",
    "SRS294940", "SRR413684",
    "SRS294941", "SRR413685",
    "SRS294942", "SRR413686",
    "SRS294943", "SRR413687",
    "SRS294944", "SRR413688",
    "SRS294945", "SRR413689",
    "SRS294946", "SRR413690",
    "SRS294947", "SRR413691",
    "SRS294948", "SRR413692",
    "SRS294949", "SRR413693",
    "SRS294950", "SRR413694",
    "SRS294951", "SRR413695",
    "SRS294952", "SRR413696",
    "SRS294953", "SRR413697",
    "SRS294954", "SRR413698",
    "SRS294955", "SRR413699",
    "SRS294956", "SRR413700",
    "SRS294957", "SRR413701",
    "SRS294958", "SRR413702",
    "SRS294959", "SRR413703",
    "SRS294960", "SRR413704",
    "SRS294961", "SRR413705",
    "SRS294962", "SRR413706",
    "SRS294963", "SRR413707",
    "SRS294964", "SRR413708",
    "SRS294965", "SRR413709",
    "SRS294966", "SRR413710",
    "SRS294967", "SRR413711",
    "SRS294968", "SRR413712",
    "SRS294969", "SRR413713",
    "SRS294970", "SRR413714",
    "SRS294971", "SRR413715",
    "SRS294972", "SRR413716",
    "SRS294973", "SRR413717",
    "SRS294974", "SRR413718",
    "SRS294975", "SRR413719",
    "SRS294976", "SRR413720",
    "SRS294977", "SRR413721",
    "SRS294978", "SRR413722",
    "SRS294979", "SRR413723",
    "SRS294980", "SRR413724",
    "SRS294981", "SRR413725",
    "SRS294982", "SRR413726",
    "SRS294983", "SRR413727",
    "SRS294984", "SRR413728",
    "SRS294985", "SRR413729",
    "SRS294986", "SRR413730",
    "SRS294987", "SRR413731",
    "SRS294988", "SRR413732",
    "SRS294989", "SRR413733",
    "SRS294990", "SRR413734",
    "SRS294991", "SRR413735",
    "SRS294992", "SRR413736",
    "SRS294993", "SRR413737",
    "SRS294994", "SRR413738",
    "SRS294995", "SRR413739",
    "SRS294996", "SRR413740",
    "SRS294997", "SRR413741",
    "SRS294998", "SRR413742",
    "SRS294999", "SRR413743",
    "SRS295000", "SRR413744",
    "SRS295001", "SRR413745",
    "SRS295002", "SRR413746",
    "SRS295003", "SRR413747",
    "SRS295004", "SRR413748",
    "SRS295005", "SRR413749",
    "SRS295006", "SRR413750",
    "SRS295007", "SRR413751",
    "SRS295008", "SRR413752",
    "SRS295009", "SRR413753",
    "SRS295010", "SRR413754",
    "SRS295011", "SRR413755",
    "SRS295012", "SRR413756",
    "SRS295013", "SRR413757",
    "SRS295014", "SRR413758",
    "SRS295015", "SRR413759",
    "SRS295016", "SRR413760",
    "SRS295017", "SRR413761",
    "SRS295018", "SRR413762",
    "SRS295019", "SRR413763",
    "SRS295020", "SRR413764",
    "SRS295021", "SRR413765",
    "SRS295022", "SRR413766",
    "SRS295023", "SRR413767",
    "SRS295024", "SRR413768",
    "SRS295025", "SRR413769",
    "SRS295026", "SRR413770",
    "SRS295027", "SRR413771",
    "SRS295028", "SRR413772",
    "SRS295029", "SRR413773"
  )
  
  combined_qin <- full_join(temp_join, Qin_key, by= "SRASample")
  
  combined_qin <- combined_qin %>%
    rename(PatientID = 'SRR') %>%
    mutate(dataset = rep("Qin", nrow(combined_qin))) %>%
    rename(type_of_disease = disease) %>%
    rename(disease = study_condition) %>%
    rename(CHOL = cholesterol) %>%
    rename(HDL = "HDL (mmol/L)") %>%
    rename(LDL =  "LDL (mmol/L)") %>%
    rename(TGL = "TG (mmol/L)")
  
  combined_qin$disease[combined_qin$disease == "IGT"] <- 1
  combined_qin$disease[combined_qin$disease == "T2D"] <- 1
  combined_qin$disease[combined_qin$disease == "control"] <- 0
  combined_qin$Ethnicity = combined_qin$country
  
  qin_curated <- combined_qin %>% select(dataset, PatientID, sampleID, disease, age, Gender, country, Ethnicity, BMI, CHOL, HDL, LDL, TGL )
  
  ## KARLSSON ## 
  GET("https://static-content.springer.com/esm/art%3A10.1038%2Fnature12198/MediaObjects/41586_2013_BFnature12198_MOESM507_ESM.xlsx", 
      write_disk(tf <- tempfile(fileext = ".xlsx")))
  more_karlsson_data <- read_excel(tf, 1L, skip = 1)
  more_karlsson_data <- more_karlsson_data %>% 
    mutate(to_merge = rep("S", nrow(more_karlsson_data))) %>%
    rename(sampleID = "Sample ID") %>%
    unite(sampleID, c("to_merge","sampleID"), sep = "") %>%
    rename(age = "Age (years)") %>%
    filter(!is.na(age)) %>%
    rename(BMI = "BMI (kg/m2)") %>%
    rename(CHOL = "Cholesterol (mmol/L)") %>%
    rename(CR = "Creatinine (µmol/L)") %>%
    rename(HDL = "HDL (mmol/L)") %>%
    rename(HSCRP = "hsCRP (mg/L)") %>%
    rename(LDL = "LDL (mmol/L)") %>%
    rename(TGL = "Triglycerides (mmol/L)")
  
  
  bugs_list <- curatedMetagenomicData("KarlssonFH_2013.metaphlan_bugs_list.stool", dryrun=FALSE)
  karlsson_data <- pData(bugs_list[[1]])
  karlsson_data <- karlsson_data %>% 
    mutate(dataset = rep("Karlsson", nrow(karlsson_data))) %>%
    rename(PatientID = NCBI_accession) %>%
    rename(sampleID = subjectID) %>% 
    rename(type_of_disease = disease) %>%
    rename(disease = study_condition) %>%
    rename(Gender = gender)
  
  karlsson_data$disease[karlsson_data$disease == "IGT"] <- 1
  karlsson_data$disease[karlsson_data$disease == "T2D"] <- 1
  karlsson_data$disease[karlsson_data$disease == "control"] <- 0
  karlsson_data$Ethnicity = karlsson_data$country
  
  karlsson_final <- full_join(more_karlsson_data, karlsson_data, by = "sampleID" )
  karlsson_final <- karlsson_final %>%
    rename(age = age.x)
  
  karlsson_curated <- karlsson_final %>% select(dataset, PatientID, sampleID, disease, age, Gender, country, Ethnicity, BMI, CHOL, CR, HDL, HSCRP, LDL, TGL )
  
  ## HMP ## 
  sample_data <- read_tsv("https://storage.googleapis.com/gbsc-gcp-project-ipop_public/HMP/clinical_tests/clinical_tests.txt")
  sample_data$VisitID <- as.character(sample_data$VisitID)
  sample_data <- sample_data %>% separate(VisitID, into = c("PatientID", "VisitNum"), sep = "-")
  
  # Patient Data  
  GET("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1236-x/MediaObjects/41586_2019_1236_MOESM3_ESM.xlsx", write_disk(tf <- tempfile(fileext = ".xlsx")))
  patient_data <- read_excel(tf, 1L)
  patient_data <- patient_data %>% rename("PatientID" = "SubjectID")
  
  #Combining Patient with Sample Data 
  combined_HMP_data <- full_join(sample_data, patient_data, by = "PatientID")
  
  # Convert Columns to correct values 
  combined_HMP_data <- combined_HMP_data %>% rename(HDL_old = HDL) %>% 
    mutate(HDL = HDL_old * 0.0555) %>%
    mutate(dataset = rep("HMP", nrow(combined_HMP_data))) %>%
    mutate(country = rep("USA", nrow(combined_HMP_data))) %>%
    rename(sampleID = VisitNum) %>%
    rename(disease = Class)%>%
    rename(age = Adj.age)
  
  combined_HMP_data$disease[combined_HMP_data$disease == "Diabetic"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Prediabetic"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Crossover"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Control"] <- 0
  
  final_HMP <- combined_HMP_data %>% select(dataset, PatientID, sampleID, disease, age, Gender, country, Ethnicity, BMI, CHOL, CR, HDL, HSCRP, LDL, TGL )
  
  
  ## Combining all three datasets together
  final_metadata <- merge(final_HMP, karlsson_curated, all = TRUE)
  final_metadata <- merge(final_metadata, qin_curated, all = TRUE)
  saveRDS(final_metadata, 'final_metadata.rds')
  
}

main()


