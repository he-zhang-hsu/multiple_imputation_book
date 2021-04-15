# line plot for survival 

require(ggplot2)

# Your example data
dat <- structure(list(ID = 1:5, eventA = c(0L, 1L, 1L, 0L, 1L), 
    eventB = c(1L, 0L, 0L, 1L, 0L), t1 = c(7, 5, 10, 4.5, 2), t2 = c(7, 5, 10, 4.5, 
    8), censored = c(0, 0, 0, 0, 1)), .Names = c("ID", "eventA", 
    "eventB", "t1", "t2", "censored"), class = "data.frame", row.names = c(NA, -5L))

# Create event variable
dat$event <- with(dat, ifelse(eventA, "A", "B"))

# Create id.ordered, which is a factor that is ordered by t2
# This will allow the plot to be ordered by increasing t2, if desired
dat$id.ordered <- factor(x = dat$ID, levels = order(dat$t2, decreasing = T))

# Use ggplot to plot data from dat object
ggplot(dat, aes(x = id.ordered)) + 
    # Plot solid line representing non-interval censored time from 0 to t1
    geom_linerange(aes(ymin = 0, ymax = t1)) + 
    # Plot line (dotted for censored time) representing time from t1 to t2
    geom_linerange(aes(ymin = t1, ymax = t2, linetype = as.factor(censored))) +  
    # Plot points representing event
    # The ifelse() function moves censored marker to middle of interval
    geom_point(aes(y = ifelse(censored, t1 + (t2 - t1) / 2, t2), shape = event), 
        size = 4) +
    # Flip coordinates
    coord_flip() + 
    # Add custom name to linetype scale, 
    # otherwise it will default to "as.factor(censored))"
    scale_linetype_manual(name = "Censoring", values = c(1, 2), 
        labels = c("Not censored", "Interval censored")) +
    # Add custom shape scale.  Change the values to get different shapes.
    scale_shape_manual(name = "Event", values = c(19, 15)) +
    # Add main title and axis labels
    # opts(title = "Patient follow-up") + xlab("Patient ID") +  ylab("Days") + 
    # I think the bw theme looks better for this graph, 
    # but leave it out if you prefer the default theme
    theme_bw()

dat

newdata=cbind.data.frame(dat$ID, dat$t1);

newdata=dat[, c("ID", "t1", "t2")];
newdata$Event=c("Failure", "Failure", "Censored", "Failure", "Failure");
newdata$EventF=c("Censored", "Censored", "Failure", "Censored", "Censored");
newdata$t1[3]=8;
newdata$t2[5]=2;
newdata$censored=c(1,1,1,1,1);
newdata$t2=c(9, 6, 10, 7.5, 3);
newdata

# plot the right-censored data;
ggplot(newdata, aes(y = ID)) + 
    # Plot solid line representing non-interval censored time from 0 to t1
    geom_linerange(aes(xmin = 0, xmax = t1)) +
 #   geom_linerange(aes(xmin = t1, xmax = t2, linetype=as.factor(censored))) +  
    geom_point(aes(x=t1, shape=Event), size=4) +
 #   geom_point(aes(x=t2, shape=EventF), size=4) +
    scale_linetype_manual(name = "Censoring", values = c(1, 2),labels = c("Not censored", "Right censored")) +
    xlab("Time") 
    
        # Flip coordinates
# plot the complete version of right-censored data;

ggplot(newdata, aes(y = ID)) + 
    # Plot solid line representing non-interval censored time from 0 to t1
    geom_linerange(aes(xmin = 0, xmax = t1, linetype=as.factor(2))) +
    geom_linerange(aes(xmin = t1, xmax = t2, linetype=as.factor(3))) +  
    geom_point(aes(x=t1, shape=Event), size=4) +
    geom_point(aes(x=t2, shape=EventF), size=4)+
    scale_linetype_manual(name = "Completeness", values = c(1, 2),labels = c("Observed", "Not Observed")) +
    xlab("Time") 
