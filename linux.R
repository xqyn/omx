# Check if the OS is Linux
if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Linux") {
  
  #--------------------------------------------------
  cat(R.version.string, '\n')
  # Library
  library(unix)
  library(BiocParallel)
  library(doParallel)
  library(parallel)
  library(pryr)
  
  # Width of display in terminal
  options(width=190)
  
  # Function --------------------------------------------
  
  # Convert bytes to gigabytes (95% of the total size)
  bytes_to_gb <- function(bytes) {
    gb <- bytes * 0.95 / 1e9    # Only use 95% of limit
    return(gb)
    }
  # Convert gigabytes to bytes (95% of the total size)
  gb_to_bytes <- function(gb) {
    bytes <- gb * 1e9 * 0.95   # Only use 95% of limit
    return(bytes)
  }
   
  # Check environment variables for assigned resources
  assigned_cores <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")  # For Slurm
  if (assigned_cores == "") {
    assigned_cores <- Sys.getenv("PBS_NP")  # For PBS
  }
  if (assigned_cores == "") {
    assigned_cores <- Sys.getenv("LSB_DJOB_NUMPROC")  # For LSF
  }
  
  # If no specific environment variable, fall back to detectCores()
  if (assigned_cores == "") {
    assigned_cores <- parallel::detectCores(logical = FALSE)  # Only physical cores
  }


  # Setting core --------------------------------------------------
  # Setting memory and core usage
  mem_bytes <- gb_to_bytes(as.numeric(gsub("[^0-9]", "", Sys.getenv("R_MAX_VSIZE"))))
  mem <- mem_bytes / 1e9
  core <- as.numeric(gsub("[^0-9]", "", Sys.getenv("R_MEMORY_LIMIT")))
  cat('mem: ', round(mem, 2), ' GB', '\n')
  cat('core: ', core, '\n')
  
  # Set resource limits
  # Set resource limits (ensure no NA values)
  if (!is.na(mem)) {
    rlimit_as(mem)
  }
  
  if (!is.na(core)) {
    rlimit_core(core)
  }
  #print(rlimit_all())

  # Register for parallel processing
  if (is.na(core) || core < 1) {
    cat('Cannot detect assigned cores: on R ood', '\n')
  } else {
    register(MulticoreParam(core))
    registerDoParallel(core) 
    register(DoparParam(), default = TRUE)
  }
  
  
  # # Memory: Checking available memory from /proc/meminfo
  # MemTotal <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE))
  #   # Convert to GB (just once)
  # MemTotal <- bytes_to_gb(MemTotal * 1024) # Convert KB to bytes then to GB
  # cat('MemTotal: ', round(MemTotal, 2), ' GB', '\n')

  # # Memory: Checking available memory from /proc/meminfo
  # MemFree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))
  #  # Convert to GB (just once)
  # MemFree <- bytes_to_gb(MemFree * 1024) # Convert KB to bytes then to GB
  # cat('MemFree: ', round(MemFree, 2), ' GB', '\n')

  #   # Memory: Checking available memory from /proc/meminfo
  # MemAvailable <- as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
  #  # Convert to GB (just once)
  # MemAvailable <- bytes_to_gb(MemAvailable * 1024) # Convert KB to bytes then to GB
  # cat('MemAvailable: ', round(MemAvailable, 2), ' GB', '\n')
} else {
  print("This script is designed to run on Linux only.")
}


check_system_resources <- function() {
  # Get number of CPU cores
  num_cores <- parallel::detectCores()
  
  # Get memory usage from the system using free -h
  memory_info <- system("free -h", intern = TRUE)
  
  # Extract relevant memory information
  memory_lines <- strsplit(memory_info, "\\s+")  # Split by whitespace
  
  # Find the line containing the memory details
  mem_line <- memory_lines[[2]]  # The second line usually contains memory info
  
  # Extract the total and used memory from the line
  total_memory <- mem_line[2]
  used_memory <- mem_line[3]
  
  # Get memory used by R using pryr::mem_used
  r_mem_used <- mem_used()
  
  # Convert R memory used from bytes to MB
  r_mem_used_mb <- r_mem_used / (1024^2)  # Convert bytes to MB
  
  # Print the results
  # cat("Number of CPU Cores: ", num_cores, "\n")
  # cat("System Total Memory: ", total_memory, "\n")
  # cat("System Used Memory: ", used_memory, "\n")
  cat("R Memory Used: ", round(r_mem_used_mb, 2), "MB\n")  # Round to 2 decimal places
}

# Call the function
check_system_resources()