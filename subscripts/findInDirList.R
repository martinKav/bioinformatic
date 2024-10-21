# source file for finding file names in a directory that match include a certain string

findInDirList  <-function(path, wantedstrings, found = FALSE){
  fileNameVec <- vector(mode = "character")
  dirList <- list.files(path)
  for (dir in dirList)
  {
    # check if found is untrue and if the folder/ file has a name we want
    if(found !=TRUE &&  grepl( paste0(wantedstrings, collapse="|"), dir ))
    {
      # set found to true
      found=TRUE
      #check if folder
      if(dir.exists(paste0(path,"/", dir)))
      {
        fileNameVec <- append(fileNameVec, findInDirList( paste0(path, "/", dir), wantedstrings, found))
      } #if not, append the filepath to the vec
      else{
        fileNameVec <- append(fileNameVec, paste0(path,"/", dir))
      }
      found=FALSE
    }
    # if found is true we call it but dont change the value
    else if (found)
    {
      # check if its a folder and recur, if not then append the file to the vec
      if(dir.exists(paste0(path,"/", dir)))
      {
        fileNameVec <- append(fileNameVec, findInDirList( paste0(path, "/", dir), wantedstrings, found))
      }else{
        fileNameVec <- append(fileNameVec, paste0(path,"/", dir))
      }
    }
    # if found is untrue and the folder is not matching, we recur anyway with a false found and append nothing
    else
    {
      if(dir.exists(paste0(path,"/", dir)))
      {
        fileNameVec <- append(fileNameVec, findInDirList( paste0(path, "/", dir), wantedstrings, found))
      }
    }
  }
  return(fileNameVec)
}
