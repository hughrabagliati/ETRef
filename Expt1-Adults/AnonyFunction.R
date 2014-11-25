# Anonymizing function, takes a list of subject names and 
# creates a set of grep commands that can be fed into the terminal
# to replace those names with anonymized replacements

anon <- function(Names){
		
		commands <- paste("grep -r -l '",unique(Names),"' . | sort | uniq | xargs perl -e \"s/",unique(Names),"/KidSubjExpt2-",1:length(unique(Names)),"/\" -pi ", sep = "")
		writeLines(commands, "AnonCommands.txt")
		}
	
