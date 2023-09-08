library(limma)
library(statmod)
library(stringi)
library(jsonlite)

get_mod_t_test_updn <- function(lfc_matrix, g2){
  design <- model.matrix(~ 0 + g2)
  num_rep <- dim(lfc_matrix)[2] / 2
  print(num_rep)
  brep <- rep(1:num_rep, each = 2) 
  dupcor <- duplicateCorrelation(lfc_matrix, design=design, block=brep)
  print(dupcor$consensus.correlation)
  
  colnames(design) <- c("group1", "group2")
  fit1 <- lmFit(lfc_matrix, design=design, block=brep, correlation=dupcor$consensus)
  cont.matrix <- makeContrasts(group1vsgroup2="group1-group2", levels=design)
  fit2 <- contrasts.fit(fit1, cont.matrix)
  fit3 <- eBayes(fit2)
  result <- topTable(fit3, coef = 1, number = Inf, confint = TRUE, sort.by = "none")[,-4]
  return(result)
}

mod_t_test <- function(input_files, output_direc, group_type='SI', plate_col='Plate') {
  files <- list.files(path=input_files, pattern="*.txt", full.names=TRUE, recursive=FALSE)
  if (!dir.exists(output_direc)) { dir.create(output_direc, recursive = TRUE) }
  for (input_direc in files){
    n_rep_str <- stri_sub(input_direc,-8,-5)
    n_rep <- 3*2
    if (n_rep_str == "2rep"){
      n_rep <- 2*2
    } else if (n_rep_str == "4rep"){
      n_rep <- 4*2
    } else if (n_rep_str == "5rep"){
      n_rep <- 5*2
    } else if (n_rep_str == "6rep"){
      n_rep <- 6*2
    } else if (n_rep_str == "7rep"){
      n_rep <- 7*2
    } else if (n_rep_str == "8rep"){
      n_rep <- 8*2
    } else if (n_rep_str == "9rep"){
      n_rep <- 9*2
    } else{
      print("rep number warning!")
      print(n_rep)
    }
    print(n_rep)
    
    plate_str <- strsplit(input_direc, "_")[[1]]
    plate_str <- plate_str[grep(plate_col, plate_str)]
    output_direc_new <- paste(output_direc, group_type, '_mod_t_test_', plate_str, '.csv', sep='')
    
    plate_all <- read.csv(file=input_direc, header=FALSE, sep='\t')
    p_matrix <- data.matrix(plate_all[2:dim(plate_all)[2]])
    g2 <- factor(c(rep("group 1", n_rep), rep("group 2", (dim(p_matrix)[2]-n_rep))))
    #p_result <- mod.t.test(p_matrix, group = g2)
    p_result <- get_mod_t_test_updn(p_matrix, g2)
    rownames(p_result) <- plate_all$V1
    write.csv(p_result, output_direc_new)
  }
}


# Read the JSON configuration file
config <- fromJSON("config.json")

# Access the parameters from the JSON object
# Set working directory
work_dir <- config$work_dir
setwd(work_dir)

# Set metatable directory
meta_dir <- config$meta_dir

# Set output directory
output_dir <- config$output_dir

# Read the metatable that guides the column names to refer to
df_meta <- read.csv(file=meta_dir, header=TRUE, sep='\t')

plate_col = unique(df_meta$plate_column)[1]

groups <- unique(df_meta$group)
for (group in groups) {
  print(paste("Group:",group))
  group_dir = paste0(output_dir, group, '/')
  mod_t_input_dir = paste0(group_dir, 'mod_t_input/')
  mod_t_output_dir = paste0(group_dir, 'mod_t_output/')
  
  final_mod_input_dir = paste0(mod_t_input_dir, group, '_LFC_var_modt_input/')
  final_mod_output_dir = paste0(mod_t_output_dir, group, '_LFC_var_modt_output/')
                            
  print("Conducting the moderated t-test...")
  mod_t_test(final_mod_input_dir, final_mod_output_dir, group_type=group, plate_col=plate_col)
}
