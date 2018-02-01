##===============================================================================================##
##                                                                                               ##
##                              LOADING REQUIRED LIBRARIES	                                     ##
##                                                                                               ##
##===============================================================================================##

if(!require(geiger)){
  install.packages("geiger",repos="http://www.stats.bris.ac.uk/R/")
  library("geiger") ## ape and gieger to read the tree
}

if(!require(ape)){
  install.packages("ape",repos="http://www.stats.bris.ac.uk/R/")
  library("ape")
}

if(!require(igraph)){
  install.packages("igraph",repos="http://www.stats.bris.ac.uk/R/")
  library("igraph") 
} # igraph for the search

if(!require(BMhyd)){ 
  install.packages("BMhyd",repos="http://www.stats.bris.ac.uk/R/")
  library("BMhyd") 
} # loaded to use 'GetAncestor' function

if(!require(phangorn)){
  install.packages("phangorn",repos="http://www.stats.bris.ac.uk/R/")
  library("phangorn")
}

if(!require(phytools)){
  install.packages("phytools",repos="http://www.stats.bris.ac.uk/R/")
  library("phytools")
}

if(!require(optparse)){
  install.packages("optparse",repos="http://www.stats.bris.ac.uk/R/")
  library("optparse")
} #Loaded for parsing

##===============================================================================================##
##                                                                                               ##
##                                        INPUT FILES                                            ##
##                                                                                               ##
##===============================================================================================##

option_list = list(
  make_option(c("-t", "--tree"), type="character", default="newick.tre", dest="tree",
              help="Rooted tree file in newick format", metavar="Newick.tre"),
  
  make_option(c("-s", "--snp_threshold"), type="integer", default="2", dest="SNPthreshold",
              help="SNP threshold for sub-grouping", metavar="SNP"),
  
  make_option(c("-r", "--relatibility_threshold"), type="integer", default="3", dest="Rthreshold",
              help="Relatibility threshold in number of internal nodes between Sub-Groups and Singletons", metavar="Internal nodes"),
  
  make_option(c("-i", "--sample_id"), type="character", default="sample_id.csv", dest="sampleID",
              help="csv file of lane id, sample id", metavar="sample_id.csv")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tree)){
  cat_help(opt_parser)
  stop("A rooted phylogenetic tree must be declared", call.=FALSE)
}

if (is.null(opt$SNPthreshold)){
  cat_help(opt_parser)
  stop("SNP threshold for Sub-Grouping much be declared", call.=FALSE)
}

if (is.null(opt$Rthreshold)){
  cat_help(opt_parser)
  stop("Relatibility threshold for Sub-Grouping much be declared", call.=FALSE)
}

tree = read.newick(opt$tree)
thresh = opt$SNP
rthreshold = opt$Rthreshold
sample_id_file = opt$sampleID
Sample.ID = read.csv(sample_id_file,sep=",",header=T, quote="", comment.char="", stringsAsFactors=FALSE)

Date <- Sys.Date() #Retrieving date for output file naming purposes

#Resolving dichotomies by deleting all branches smaller than 0.5 i.e. zero SNPs and collapses the corresponding dichotomies into a multichotomy.
tree = di2multi(tree, 0.5)
tree_name = paste("NewickTreefile_",Date,".tre")
write.tree(tree, file=tree_name) #output the new tree

#Display Number of Taxa in Tree 
cat("\t",length(tree$tip.label)," taxa found from phylogenetic tree\n",sep="")

#Renaming the node labels; at a sequence from 1 to (number of nodes in tree), at increments of 1
tree$node.label = paste("Node",seq(1,tree$Nnode,1),sep="_")  


##===============================================================================================##
##                                                                                               ##
##                                    FUNCTION DECLARATIONS                                      ##
##                                                                                               ##
##===============================================================================================##

##===============================================================================================##
##                                                                                               ##
##                   Function: Getting distance between two nodes + distance                     ##
##                                                                                               ##
##===============================================================================================##


#Given
# graph_tree: a rooted tree of newick format (usually the "input_tree") converted into graph form by igraph
# root: a node in the graph you want to begin the path
# node: target node where you want to the path to run to
# input_tree: a rooted tree that has been parsed into this script - class "phylo"
#Return the distance between node and tip

get.path.distance <- function(graph_tree, root, node,input_tree){
  path <- shortest_paths(graph_tree, 
                         root, 
                         node, 
                         mode = "out", 
                         weights = NULL, 
                         output = "both") #Determine path between nodes
  
  #retrieve edgelist of tree in igraph form (different from the original tree)
  igraph_edgelist <- as_edgelist(graph_tree) 
  
  #define a matrix to store edges
  edge_matrix <-  matrix(ncol=2, nrow=length(path$epath[[1]])) 
  
  # retrieve edges of path and store them
  for(i in 1:length(path$epath[[1]])){ 
    
    #stating the edge; retrieved as "row" number in edgelist vector
    edge <- as.integer(path$epath[[1]][i]) 
    
    #retrieving the edge from list
    edge_matrix[i,] <- igraph_edgelist[edge,] 
  }
  #comparing two rows; if both elements match i.e. all() will call TRUE
  step1 <- Vectorize(function(x,y) { 
    all(input_tree$edge[x,] == edge_matrix[y,]) 
  })
  
  #Generate a matrix of a pairwise comparison between two matrices and their edges
  #Rows correspond to tree$edges vs cols correspond to edges in path 
  result <- outer(1:nrow(input_tree$edge), 
                  1:nrow(edge_matrix), 
                  FUN=step1) 
  #identify which row each edge is found in input_tree$edge
  edge_row <- which(apply(result, 1, any)) 
  
  path_length = 0
  
  for(i in edge_row){ #for each edge row, return its corresponding edge length and add to the path length
    edge_length = as.integer(input_tree$edge.length[i])
    path_length = path_length + edge_length
  }
  
  return(path_length)
  
} #End of Function

get.path <- function(graph_tree, root, node) #inputs: igraph.tree, root, target node,input_tree is the newick tree
{
  path <- shortest_paths(graph_tree, root, node, mode = "vpath", weights = NULL, output = "both") #Determine path between nodes
  return(path$vpath)
}

##===============================================================================================##
##                                                                                               ##
##               Function: Identify tips within zero distance of internal node                   ##
##                                                                                               ##
##===============================================================================================##


#required - subtree; igraph.tree, #graph;node, #node that is the root of the subtree ; tip,#tip define distance from said node; tree
tips.zero.distance.node <- function(node,subtree,igraph.tree,tree){
  tips_under_node <- subtree[which(subtree <= ntips)] #retrieve tips under said node
  tips_zero_distance_node <- rep(0,length(tips_under_node)) #Create empty vector for storage
  interation_0dist_vector = 0
  for(x in tips_under_node){ #distances of tips from node
    tip <- as.numeric(as.character(x))
    node_tip_dist <- get.path.distance(igraph.tree, #graph
                                       node, #node that is the root of the subtree
                                       tip,#tip define distance from said node
                                       tree) #retrieving the distances from node to tip
    #going to store the tips that have a zero into a vector	  
    if(node_tip_dist == 0){
      interation_0dist_vector <- interation_0dist_vector+1
      tips_zero_distance_node[interation_0dist_vector] <- tip
      
    }
  }
  tips_zero_distance_node <- tips_zero_distance_node[tips_zero_distance_node != 0] #clearing elements in the vector that are zero
  return(tips_zero_distance_node)
}

##===============================================================================================##
##                                                                                               ##
##             Function: Identify  distance of all times under a internal node                   ##
##                                                                                               ##
##===============================================================================================##

tips.distance.node <- function(node,subtree,igraph.tree,tree){
  
  #retrieve tips under said node
  tips_under_node <- subtree[which(subtree <= ntips)] 
  
  #Create empty vector for storage
  tips_distance_node <- rep(0,length(tips_under_node)) 
  interation_dist_vector = 0
  
  for(x in tips_under_node){ #distances of tips from node
    tip <- as.numeric(as.character(x))
    node_tip_dist <- get.path.distance(igraph.tree, #graph
                                       node, #node that is the root of the subtree
                                       tip,#tip define distance from said node
                                       tree) #retrieving the distances from node to tip
    
    #going to store the tip distances into a vector	  
    interation_dist_vector <- interation_dist_vector+1
    tips_distance_node[interation_dist_vector] <- node_tip_dist
    
  }
  return(max(tips_distance_node)) #Return the max distance 
}


###################################################################################################
##                                                                                               ##
##                                Main Subgrouping Algorithms                                    ##
##                                                                                               ##
###################################################################################################



##===============================================================================================##
##                                                                                               ##
##                               Define variables for Sub-grouping                               ##
##                                                                                               ##
##===============================================================================================##

#Number of tips
ntips<-Ntip(tree) 

# setting a Sub-Group number - will tick up as more clusters are identified
sgnum <- 0 

# cluster assignment - rep replicates the values found in x. rep(x, ...); generating a vector the size of ntips+tree$Nnode corresponding to the number of positions within the tree data
assign <- rep(0,ntips+tree$Nnode) 

# convert tree in igraph form - Takes a 2 column matrix edge list, where each row is an existing edge 
igraph.tree <- graph.edgelist(tree$edge) 

#Depth First Search traverses a graph, begins at a "root" vertex and tries to go quickly as far from as possible 
#Order - return the DFS ordering of the vertices; dist - return the distance from the root of the search tree
dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',order=TRUE,dist=TRUE) 

##===============================================================================================##
##                                                                                               ##
##                                  Investigating edge distances                                 ##
##                                                                                               ##
##===============================================================================================##

#Traverse the tree in depth first order - starting at the root
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    
    # Skip leaves/tips
    #If the number of the "node" is less than ntips+1 it is actually a leaf/tip 
    #i.e. phytools stores tips of trees as 1:n (n = no. of tips) then nodes as n+1:n+m (m = no. of nodes)
    if(node < ntips+1) {
	    cat(paste(tree$tip.label[node]," skipped as this is a leaf."), sep="\n")
 	    next
    } else
      
    {
      #if the number of the "node" is equal or more than ntips+1 it is actually a node; proceed with below
      cat("", paste(tree$node.label[node-ntips]," will be investigated."), sep="\n") 
    }
    
    #Assign a subgroup number containing the said node & all the tips connected in the 'out' direction from the root 
    if(assign[node]<=0){ 
      
      #Said node becomes the root; dfs will search out and notes the order of "nodes" (which includes internal nodes + tips); display only the $order and store it
      subtree <- graph.dfs(igraph.tree,node,neimode='out',unreachable=FALSE)$order 
      
      #Checks within the igraph.vs called subtree which values do not equal to NA and stores them under subtree i.e. remove all NAs
      subtree <- subtree[!is.na(subtree)] 
      
      #Define maximum distance amongst all tips from said
      tips_dist_node <- tips.distance.node(node,subtree,igraph.tree,tree)	
  
  		#Define the tips connected to node at zero distance
  		tips_zero_distance_node <- tips.zero.distance.node(node,subtree,igraph.tree,tree)
  		
  		#Logical sub-grouping numbering; composed of  two parameters - threshold && zero distance
  		
  		#Logical: tip distances from node below or equal to threshold?
      if(tips_dist_node<=thresh){ 
        #Tick up the clustering number
        sgnum <- sgnum+1 
        
        #Record/assign the cluster number to the "assign" vector in the positions specified by the "subtree" variable
        assign[subtree] <- sgnum 
        
        #Display that a Sub-Group has been identified and assigned
        cat(paste("Threhold met - Sub-Group Number assigned: ",sgnum), sep="\n")
      }
  		
  		#If said internal node has downstream tips over the declared SNP threshold assess tips directly connected to said node
  		#Assign Sub-group if there are tips are zero distance to node; Must not be a singleton
  		else if(sum(tips_zero_distance_node)>0 && tips_zero_distance_node[1] != sum(tips_zero_distance_node)){ 
        
  		  sgnum <- sgnum+1
        
        #assign sgnum to tips at node in zero distance 
        assign[tips_zero_distance_node] <- sgnum
        
        #assign sgnum to the node
        assign[node] <- sgnum 
        
        #Display that a Sub-Group has been identified and assigned
        cat(paste("Internal sub-group detected - Sub-group number assigned: ",sgnum), sep="\n")
        
  		}
  		
    } #End of If statement
    
} # End of for loop



##===============================================================================================##
##                                                                                               ##
##                            Identification of Major Sub-Grouping                               ##
##                                                                                               ##
##===============================================================================================##
  
iteration = 0
sg_intersect_list <- list()

for(i in 1:sgnum){
  #which elements are found for a sgnum
  sgnum_elements <- which(assign==i)
  sgnum_elements_tips <- sgnum_elements[which(sgnum_elements< ntips+1)]
  iteration = 0
  sg_element_ancestors_list <- list()
  
  #Retrieve the internal nodes of each member of said subgroup
  for(x in sgnum_elements_tips){
    
    #generate a list of each SG element
    sg_element_ancestors <- Ancestors(tree, x, type=c("all")) 
    iteration = iteration + 1
    
    #Store the list of ancestors
    sg_element_ancestors_list[iteration] <- list(sg_element_ancestors) 
  }
  
  #Collapse the list into intersects and store for that sgnum
  sub_intersect <- Reduce(intersect,sg_element_ancestors_list)
  
  #Store the list of intersecting ancestors
  sg_intersect_list[i] <- list(sub_intersect) 
}

#combinations in columns
combinations <- combn(length(sg_intersect_list), 2) 

#Take the ancestor nodes of each Sub-Group and setup a paired-wise comparison  - list of lists of the combinations
ll <- combn(sg_intersect_list, 2 , simplify = FALSE ) 

 # Pair-wise Comparison Function! Intersect the list elements - find the intersect &  the length  
 out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )

 #which Pair-wise comparisons have more than two internal nodes
 
 #Declare an empty variable to store wanted comparisons
 major_subgroup_intersect <- list() 
 
 #which lists in the Pair-wise comparison list of lists have equal or more than internal node relatibility threshold
 major_subgroup_comb <- which(out >= rthreshold) 
 
 #Extracting the subgroup number for the above
 for(i in 1:length(major_subgroup_comb)){
   combin_element <- combinations[,major_subgroup_comb[i]]
   major_subgroup_intersect[[i]] <- combin_element
 }

 #Combination Function! need to combine lists with overlapping sgnum; Goal: If there are overlapping numbers in differents lists combine them
 for(i in seq_along(major_subgroup_intersect)[-length(major_subgroup_intersect)]) {
   if(length(intersect(major_subgroup_intersect[[i]], major_subgroup_intersect[[i+1]])) > 0) {
     
     major_subgroup_intersect[[i+1]] <- sort.int(unique(c(major_subgroup_intersect[[i]],
                                                          major_subgroup_intersect[[i+1]])))
     major_subgroup_intersect[[i]] <- as.list(NULL)
     
   }
 }
 
#Take only the rows with something in it - The major subgroup list is then created
major_subgroup <- Filter(function(x) length(x) > 0, major_subgroup_intersect)

#Number of Major SubGroup
numbMajorSubgroup <- length(major_subgroup) 

#Retrieving the SubGroup numbers for each Major SubGroup, retrieve the SubGroup's tips and storing them in a vector according to the tree
maj_subgroup_mems = NULL
maj_subgroup_assign <- rep(0,ntips+tree$Nnode)

for(x in 1:length(major_subgroup)){
 maj_subgroup_mems = NULL
 
 #Within the Major Sub-Group, call the tips of each Sub-Group and store them in "maj_subgroup_mems"
 for(i in major_subgroup[[x]]){ 
   maj_subgroup_mems <- c(maj_subgroup_mems,which(assign==i))
 }
 
 #Similar to the "assign" vector, for each tip noted with the major cluster in the positions of the tips
 maj_subgroup_assign[maj_subgroup_mems] <- x 

}

## recording the SGnums and their respective major Sub-Group for future reference
sgnum.major.subgroup.list <- matrix(ncol=2,nrow=length(unlist(major_subgroup))) 
sgnum.major.subgroup.list.interation = 0

#each list is a major cluster i.e. [[1]] is the 1st major cluster
for(i in 1:length(major_subgroup)){
  
  #extracting each SGnum 
  for(j in major_subgroup[[i]]){
    
    sgnum.major.subgroup.list.interation = sgnum.major.subgroup.list.interation +1 
    
    #appending each row with sgnum and respective major cluster
    sgnum.major.subgroup.list[sgnum.major.subgroup.list.interation,] <-c(j,i) 
    
  }
  
}

##===============================================================================================##
##                                                                                               ##
##                 Identification of Major Sub-Grouping: Singleton Analysis                      ##
##                                                                                               ##
##===============================================================================================##   


# Check all singletons if they are related to isolates in a Sub-Group; if that cluster is apart of the Major Sub-Group, assign singleton to Major Sub-Group

#Step 1: Retrieve ancestor nodes of each singleton and store them
singleton_ancestor_list <- NULL

#Singletons are those without an assign sgnum
singleton_elements <- which(assign[1:ntips]==0) 

#Declare how may singletons 
nsingleton <- length(singleton_elements) 

#combinations in columns
combinations <- combn(length(sg_intersect_list), 2) 

#Retrieve ancestors for said singleton & store into a list
for(i in 1:nsingleton){ 
 said_singleton = singleton_elements[i]
 
 #generate a list of each sgnum taking the 1st element
 singleton_ancestors <- Ancestors(tree, said_singleton, type=c("all")) 
 
 #list in order of singleton number
 singleton_ancestor_list[i] <- list(singleton_ancestors) 
 
}

#Step 2: Create a basis to compare the ancestor nodes of each singleton with those of each sgnum

#Pairwise comparison combinations are based per singleton against all sgnums; 
#ancestor nodes of said singleton will be listed in var2, in combination with the number sgnums and their respective ancestor nodes
#i.e each singleton comparison is in a set of rows the size the number of sgnums identified
#e.g. first singleton is row 1 to row(no of cums), second singleton is the next set of rows after to the number of sgnums
result.df <- expand.grid(sg_intersect_list,singleton_ancestor_list) 

#for each singleton retrieve the "block" of said singleton + sgnum combinations 
for(i in 1:nsingleton){
  
  #Define the block
  singleton.sgnum.length <-rep(0,sgnum) 
  single.block.start = (sgnum * i) - sgnum +1
  single.block.end = i * sgnum
  single.block = result.df[single.block.start:single.block.end,]
  
  #x is the sgnum number
  for(x in 1:nrow(single.block)){
    
    #Determine the intersect for a said row i.e. Singleton vs sgnum(which is row number)
    int <- Reduce(intersect,lapply(single.block,"[[",x)) 
    
    #Storing how many internal nodes of a Sub-Group intersect with said singleton over or equal to set Relibility Threshold
    if(length(int)>=rthreshold){ 
      
      singleton.sgnum.length[[x]] <- length(int) 
    
    }
    
    #which sgnum overlap  equal or more in internal nodes with said singleton no
    sgnum.for.singleton.int <- which(singleton.sgnum.length>=rthreshold) 
    
    #singleton needs to have intersected with a sgnum && the sgnum needs to be have a major cluster to call it
    if(length(sgnum.for.singleton.int)>0 && any(sgnum.major.subgroup.list[,1] %in% sgnum.for.singleton.int) ){
     
      #given the list of each sgnum and their corresponding major Sub-Group - retrieve corresponding major Sub-group
      majSB.of.singleton <- sgnum.major.subgroup.list[which(sgnum.major.subgroup.list[,1] %in% sgnum.for.singleton.int),#rows which desired sgnums
                                                      2 #2nd column is where the major cluster is stored
                                                      ]

      #If a singleton has a unique major cluster, store it in the major cluster assigning vector
      #results are a vector calling a Major cluster multiple times i.e. singleton can have same ancestral nodes with multiple sgnums
      
      #check if major cluster is unique, should have only one calling of a major cluster number.
      if(length(unique(majSB.of.singleton)) ==1){
        
        #assign the major cluster to the singletons in maj_subgroup_assign
        majclust <- unique(majSB.of.singleton)
        
        #convert singleton number into the actual singleton tip number
        singleton.for.assign <- singleton_elements[i] 
        
        #assign the major cluster identified to the singleton
        maj_subgroup_assign[singleton.for.assign] <- majclust
        
      }else{
        cat(paste("Error: Two Major Subgroups identified for a singleton ",unique(majSB.of.singleton)), sep="\n")
      }
      
    } #End of If Statement: Major Cluster of a Singleton Identification
    
  }
  
}

 ## collating variables from above
 ans <- list(subgroupmems=assign, #Assigned subgroup number for each leaf/tip
             allcsize=table(assign), #sizes for each subgroup; including both leaves and internal nodes
             leafclustsize=table(assign[1:ntips]), #sizes for each subgroup; leaves/tips only
             ntips=ntips, #stating the number of tips in the tree
             threshold=thresh, #stating the maximum distance that will define a claded
             majorsubgroupmems=maj_subgroup_assign) 
 ans$subgroupmems <- ans$subgroupmems[1:ntips] #subsetting the "membership" to only include tips/leaves
 ans$majorsubgroupmems <- ans$majorsubgroupmems[1:ntips]
 
 #Number of singletons after sub-grouping
 remaining_singletons <- which(assign[1:ntips] == 0)


##===============================================================================================##
#                                                                                                 #
#                                    Saving subgrouping results                                   #
#                                                                                                 #
##===============================================================================================##

  
## Compiling all meta data, subgrouping and major sub-groups for output.

if (exists(opt$sampleID)){
  majorsubgroupTable = mat.or.vec(length(tree$tip.label),4)
  colnames(majorsubgroupTable)=c("Taxa","Sample.ID","Sub-group","Major.Sub-group") #Giving column names
  majorsubgroupTable[,1] = tree$tip.label #Assign the 1st column of said matrix or vector to the name of the tips from the tree

  filter_lanes <- match(majorsubgroupTable[,1],Sample.ID[,1]) #matches in majorcluster 
  Sample.ID.order <- Sample.ID[filter_lanes, ] #reordering sample id via tree order
  majorsubgroupTable[,2] = Sample.ID.order[,2] #Assign the 2nd column of said matrix or vector to cluster it has been assigned
  
  majorsubgroupTable[,3] = ans$subgroupmems  #Assign the 2nd column of said matrix or vector to cluster it has been assigned
  majorsubgroupTable[,4] = ans$majorsubgroupmems #Assign the 3rd column of said matrix or vector to Major cluster it has been assigned
  majorzero_elems = which(majorsubgroupTable[,3]==0) #store positions/elements in the vector are equal to zero i.e. were not assigned a cluster number
  majorsubgroupTable[majorzero_elems,3]="singleton_" #Said positions in vector are replaced with "singleton" 
  singletonsii = which(majorsubgroupTable[,3]=="singleton_") #store positions/elements in the vector have singleton i.e. were not assigned a cluster number
  majorsubgroupTable[singletonsii,3]=paste("singleton_",seq(1,length(singletonsii),1),sep="") #rename each position stated in "ii" with singleton plus a unique number using seq(1 to length of ii, by increments of 1)
  
} else {
  majorsubgroupTable = mat.or.vec(length(tree$tip.label),3)
  colnames(majorsubgroupTable)=c("Taxa","Sub-group","Major.Sub-group") #Giving column names
  majorsubgroupTable[,1] = tree$tip.label #Assign the 1st column of said matrix or vector to the name of the tips from the tree
  majorsubgroupTable[,2] = ans$subgroupmems #Assign the 2nd column of said matrix or vector to cluster it has been assigned
  majorsubgroupTable[,3] = ans$majorsubgroupmems #Assign the 3rd column of said matrix or vector to Major cluster it has been assigned
  majorzero_elems = which(majorsubgroupTable[,2]==0) #store positions/elements in the vector are equal to zero i.e. were not assigned a cluster number
  majorsubgroupTable[majorzero_elems,2]="singleton_" #Said positions in vector are replaced with "singleton" 
  singletonsii = which(majorsubgroupTable[,2]=="singleton_") #store positions/elements in the vector have singleton i.e. were not assigned a cluster number
  majorsubgroupTable[singletonsii,2]=paste("singleton_",seq(1,length(singletonsii),1),sep="") #rename each position stated in "ii" with singleton plus a unique number using seq(1 to length of ii, by increments of 1)
}


majSubgroupOutput = paste(thresh,"_SNPs_SB_threshold_",numbMajorSubgroup,"_MajSubgroups__",Date,".txt",sep="")

write.table(majorsubgroupTable,file=majSubgroupOutput,sep=",",col.names=T,row.names=F,quote=F)

##===============================================================================================##
#                                                                                                 #
#                              Display Sub-Grouping stats to terminal                             #
#                                                                                                 #
##===============================================================================================##

cat("",paste("Number of Sub-groups identified: ",sgnum),
    paste("Number of Major Subgroups identified: ",length(major_subgroup)),
    paste("Number of Singletons remain: ",length(remaining_singletons)),
    sep = "\n")

#Display number of samples in major clusters and cluster in each major cluster
for(i in 1:numbMajorSubgroup){
  numbMajorSubgroupMems = length(which(majorsubgroupTable[,3]==i))
  numbMajorSubgroupsgnum = length(which(sgnum.major.subgroup.list[,2]==i))
  cat(paste("Major Sub-Group ",i,"is composed of ",numbMajorSubgroupsgnum," Sub-Groups & ",numbMajorSubgroupMems," isolates."),sep="\n")
}

