#LHM
#11/27/18
#ATiB Module 7 Code v1

#A few comments, and some differences between this model and the model used in the paper:

#The model attempts to simplify interactions/calculations somewhat, rather than using all of the formulae present in the paper.

#The model populates a 200 x 200 grid (like that in the paper) with a selection of bacteria from three strains.
#In this case, all cells in the grid are populated; the model could be modified to start with some unoccupied (0) cells,
#as desired.

#Strain 1 kills strain 2 kills strain 3 kills strain 1. Strain 1 secretes a factor that leads to the decay of the
#product emitted by strain 2 (that kills strain 3), and so on.

#During each cycle, a "kill radius" and a "decay radius" are set. The kill radius represents how far the effects
#of a deadly bacteria can travel (e.g., if the kill radius is 3, a colony from strain 1 can contribute to the death of
#any strain 2 colony found within that 3-square radius). The "decay radius" represents the range of suppression (e.g.,
#strain 1 contributes to stopping strain 2 from killing strain 3 if strain 1 is found within the decay radius distance).
#The decay minimum and kill minimum represent how many colonies must lie within a specified radius, at minimum, to
#trigger the killing or decaying effect.

#At the end of each round, spots with dead bacteria can be filled in according to a weighted formula that depends on
#the surrounding colonies and the growth rate for each strain.

#Populations of each strain at the start of each cycle are graphed.



#Required packages:
library(ggplot2)

#Creating a random starting grid for all three strains. Filling every spot with these values.
#Alternatively, could adjust this to include a chance of some zeroes so that not every spot contains a value.
start_grid <- matrix(sample.int(3, 40000, replace = TRUE), 200, 200) #Adjust number of values generated and dimensions as desired

kill_radius <- 3
decay_radius <- 3
kill_min <- 4 #Number of bugs in range needed to kill, in absence of neutralization
decay_min <- 16 #Number of bugs in range needed to neutralize

#Desired cycle information:
desired_cycles <- 50

#Setting up the data frame to use for plotting later
counts_by_cycle <- matrix(nrow = desired_cycles, ncol = 4, dimnames = list(c(1:desired_cycles), c("Cycle", "Strain_1", "Strain_2", "Strain_3")))
counts_df <- as.data.frame(counts_by_cycle)


for (c in 1:desired_cycles) {
  
  #Update the table I'm planning on using to plot stuff at the beginning of each cycle, then do the operations
  counts_df$Cycle[c[1]] <- c
  counts_df$Strain_1[c[1]] <- sum(start_grid == 1)
  counts_df$Strain_2[c[1]] <- sum(start_grid == 2)
  counts_df$Strain_3[c[1]] <- sum(start_grid == 3)
  
  
  
  #Make a duplicate of the grid before making changes. Refer back to the starting grid, then make changes in the duplicated grid. At the end, once every row and column has been looked at, can update the start_grid and begin again.
  grid_dup <- start_grid
  
  for (i in 1:dim(start_grid)[1]) {
    for (j in 1:dim(start_grid)[2]) {
      
      cell_id <- start_grid[i,j] #Keep track of the identity of a particular cell in the grid
      
      #Modifying the kill boundaries to stop at the edges as a starting point. This approach does not wrap effects as seen in the paper.
      #Could adjust in a later version to do wrapping:
      
      i_lower_kill <- (i - kill_radius) #Set the lower bound, then adjust it based on the conditions listed below, as needed
      
      if ((i-kill_radius)<1) {
        i_lower_kill <- 1
      }
      
      j_lower_kill <- (j - kill_radius) #Set the lower bound, then adjust it based on the conditions listed below, as needed
      
      if ((j-kill_radius)<1) {
        j_lower_kill <- 1
      }
      
      i_upper_kill <- (i + kill_radius) #Set the upper bound, then adjust it based on the conditions listed below, as needed
      if ((i+kill_radius)>dim(start_grid)[1]) {
        i_upper_kill <- dim(start_grid)[1]
      }
      
      j_upper_kill <- (j + kill_radius) #Set the upper bound, then adjust it based on the conditions listed below, as needed
      if ((j+kill_radius)>dim(start_grid)[2]) {
        j_upper_kill <- dim(start_grid)[2]
      }
      
      local_kill_view <- start_grid[(i_lower_kill):(i_upper_kill), (j_lower_kill):(j_upper_kill)]
      
      kill_1_counts <- sum(local_kill_view == 1)
      kill_2_counts <- sum(local_kill_view == 2)
      kill_3_counts <- sum(local_kill_view == 3)
      
      #Case where 1 kills 2 kills 3 kills 1 and 1 inhibits 2 from killing 3 and 2 inhibits 3 from killing 1 and 3 inhibits 1 from killing 2. For decay, want to take decay into account if there is a killer close enough to be an issue (otherwise, leave the position the same)
      
      #If there is enough killer around to potentially kill, check whether there is a decay factor to offset it
      if ((cell_id == 1 && kill_3_counts >= kill_min) || (cell_id == 2 && kill_1_counts >= kill_min) || (cell_id == 3 && kill_3_counts >= kill_min)) {
        
        
        #Modifying the decay boundaries to stop at the edges as a starting point. This approach does not wrap as seen in the paper.
        #Could adjust in a later version to do wrapping:
        
        i_lower_decay <- (i - decay_radius) #Set the lower bound, then adjust it based on the conditions listed below, as needed
        
        if ((i-decay_radius)<1) {
          i_lower_decay <- 1
        }
        
        j_decay_kill <- (j - decay_radius) #Set the lower bound, then adjust it based on the conditions listed below, as needed
        
        if ((j-decay_radius)<1) {
          j_lower_decay <- 1
        }
        
        i_upper_decay <- (i + decay_radius) #Set the upper bound, then adjust it based on the conditions listed below, as needed
        if ((i+decay_radius)>dim(start_grid)[1]) {
          i_upper_decay <- dim(start_grid)[1]
        }
        
        j_upper_decay <- (j + decay_radius) #Set the upper bound, then adjust it based on the conditions listed below, as needed
        if ((j+decay_radius)>dim(start_grid)[2]) {
          j_upper_decay <- dim(start_grid)[2]
        }
        
        local_decay_view <- start_grid[(i_lower_decay):(i_upper_decay), (j_lower_decay):(j_upper_decay)]
        
        decay_1_counts <- sum(local_decay_view == 1)
        decay_2_counts <- sum(local_decay_view == 2)
        decay_3_counts <- sum(local_decay_view == 3)
        
        #If there is enough of the decay factor in the area, keep the cell the same
        if ((cell_id == 1 && decay_2_counts >= decay_min) || (cell_id == 2 && decay_3_counts >= decay_min) || (cell_id == 3 && decay_1_counts >= decay_min)) {
          grid_dup[i,j] <- cell_id
        }
        
        #If there isn't enough of the decay factor in the area, cell gets killed and is replaced with zero, to be changed at the end of the round
        else {
          grid_dup[i,j] <- 0
        }
      }
      
      #If there isn't enough killer in the vicinity, leave the bacteria strain at that position the same
      else {
        grid_dup[i,j] <- cell_id #Value stays the same if nothing kills it. Calculate updated changes at the end.
      }
      
    }
  }
  
  grid_fill <- grid_dup #Making another copy of the grid to fill in
  
  #Growth rates for the different strains, used for probabilities to fill in each spot. Can adjust these as desired.
  growth_1 <- 1
  growth_2 <- 1
  growth_3 <- 1
  
  #Go through each value and check if it's zero. If it is, look at the surrounding cells. Pick a value to fill in
  #based on the growth rates and the relative number of each pixel in the surrounding region. Slightly different
  #from the paper's formula, but seems to work well for growth in a full bacterial environment.
  
  fill_radius <- 1 #Can look at cells around each empty bacterial location. Default this to 1 to just look at the directly adjacent colonies.
  
  for (i in 1:dim(grid_dup)[1]) {
    for (j in 1:dim(grid_dup)[2]) {
      if (grid_dup[i,j] == 0) {
        
        i_lower_fill <- (i - fill_radius)
        
        if ((i-fill_radius)<1) {
          i_lower_fill <- 1
        }
        
        j_lower_fill <- (j - fill_radius)
        
        if ((j-fill_radius)<1) {
          j_lower_fill <- 1
        }
        
        i_upper_fill <- (i + fill_radius) #Set the upper bound, then adjust it as needed
        if ((i+fill_radius)>dim(grid_dup)[1]) {
          i_upper_fill <- dim(grid_dup)[1]
        }
        
        j_upper_fill <- (j + fill_radius) #Set the upper bound, then adjust it as needed
        if ((j+fill_radius)>dim(grid_dup)[2]) {
          j_upper_fill <- dim(grid_dup)[2]
        }
        
        local_fill_view <- start_grid[(i_lower_fill):(i_upper_fill), (j_lower_fill):(j_upper_fill)]
        
        #This will get surrounding counts. Center is zero.
        fill_1_counts <- sum(local_fill_view == 1)
        fill_2_counts <- sum(local_fill_view == 2)
        fill_3_counts <- sum(local_fill_view == 3)
        
        fill_1_weighted <- fill_1_counts * growth_1
        fill_2_weighted <- fill_2_counts * growth_2
        fill_3_weighted <- fill_3_counts * growth_3
        
        grid_fill[i,j] <- sample(3, 1, replace = TRUE, prob = c(fill_1_weighted, fill_2_weighted, fill_3_weighted))
        
      }
    }
  }
  
  start_grid <- grid_fill #Update the starting grid before the next cycle begins
  
}

#Plot results
ggplot(counts_df, aes(Cycle)) + 
  geom_line(aes(y = Strain_1, colour = "Strain_1")) + 
  geom_line(aes(y = Strain_2, colour = "Strain_2")) +
  geom_line(aes(y = Strain_3, colour = "Strain_3")) +
  xlab("Cycle") +
  ylab("Population Per Strain") +
  ggtitle("Simulated Strain Population Over Time") +
  labs(colour='Strain')
