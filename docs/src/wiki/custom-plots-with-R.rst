Demix output manipulation and plotting using R
-------------------------------------------------------------------------------

Here we discuss freyja ``demix`` output manipulation and plotting using R programming
language.

1. Change freyja ``demix`` output to a dataframe format.
    Code credit for parsing the demix output: https://github.com/a-roguet
.. code:: R

    library(data.table)
    library(tidyverse)
    library(lubridate)

    results<-read.table("test_sweep.tsv", fill = TRUE, sep = "\t", h=T)
    results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
    results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

    # code credit: https://github.com/a-roguet
    # Summarized data
    summarized<-as.data.frame(setDT(tstrsplit(as.character(results$summarized), " ", fixed=TRUE))[]) # Extract Summarized data
    summarized$sample<-results$X

    for(i in 1:((ncol(summarized)-1)/2)){
      if(i==1){
        summarized.final<-summarized[,c(ncol(summarized),1:2)]
      } else {
        start=i*2-1; end=i*2
        summarized.final<-rbind(summarized.final, setNames(summarized[,c(ncol(summarized), start:end)], names(summarized.final)))
      }
    }
    summarized.final<-summarized.final[complete.cases(summarized.final), ]
    names(summarized.final)<-c("Sample", "lineage", "abundance")


    # Sublineages data
    for(i in 1:nrow(results)){
      lineages.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 3]), " ", fixed=TRUE))[]))
      abundances.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 4]), " ", fixed=TRUE))[]))
      sample.temp<-rep(results[i, 1], nrow(lineages.temp))
      if(i==1){
        sublineages.final<-cbind(sample.temp, lineages.temp, abundances.temp)
      } else {
        sublineages.final<-rbind(sublineages.final, cbind(sample.temp, lineages.temp, abundances.temp))
      }
    }
    names(sublineages.final)<-c("Sample", "sublineage", "abundance")


2. Read in time metadata information and merge it with abundance dataframe

.. code:: R

    time_metadata <- fread("sweep_metadata.csv")
    combined_lineage_abundance_time_data <- time_metadata %>%
      inner_join(summarized.final, by = "Sample")
    combined_sublineage_abundance_time_data <- time_metadata %>%
      inner_join(sublineages.final, by = "Sample")
3. Create lineage prevalence stacked bar plots grouped by month interval

.. code:: R

     combined_lineage_abundance_time_data %>%
      mutate(month = month(mdy(sample_collection_datetime)))%>%
      mutate(abundance = as.numeric(abundance))%>%
      group_by(lineage,month)%>%
      summarise(mean_monthly_abundance = mean(abundance)) %>%
       ggplot(grouped_interval_lineage, mapping = aes(fill= lineage, y=mean_monthly_abundance, x=month)) +
       geom_bar(position="fill", stat="identity") + theme_minimal() +ylab("Variant Prevalence") +
       ggtitle("Lineage prevalence")

4. Create sub-lineage prevalence stacked bar plots grouped by month interval

.. code:: R

      combined_sublineage_abundance_time_data %>%
      mutate(month = month(mdy(sample_collection_datetime)))%>%
      mutate(month = as.factor(month)) %>%
      mutate(abundance = as.numeric(abundance))%>%
      group_by(sublineage,month)%>%
      summarise(mean_monthly_abundance = mean(abundance))%>%
      ggplot(grouped_interval_sublineage, mapping = aes(fill=sublineage, y=mean_monthly_abundance, x=month)) +
      geom_bar(position="fill", stat="identity") + theme_minimal() +ylab("Variant Prevalence") +
      ggtitle("Sublineage prevalence") +
      theme(legend.title = element_text( size=2), legend.text=element_text(size=5))

5. Create lineage prevalence per sample plots

.. code:: R

    combined_lineage_abundance_time_data %>%
      mutate(abundance = as.numeric(abundance))%>%
      group_by(lineage, Sample)%>%
      summarise(mean_sample_abundance = mean(abundance)) %>%
      ggplot(grouped_interval_lineage, mapping = aes(fill= lineage, y=mean_sample_abundance, x=Sample)) +
      geom_bar(position="fill", stat="identity") + theme_minimal() +ylab("Variant Prevalence") +
      ggtitle("Lineage prevalence") + theme(axis.text.x = element_text(angle = 45))

6. Create stacked area plot showing lineage prevalence based on dates

.. code:: R

    combined_lineage_abundance_time_data %>%
      mutate(abundance = as.numeric(abundance))%>%
      ggplot(aes(x=sample_collection_datetime,y=abundance,group=lineage,fill=lineage)) +
      geom_area(position="fill")+ theme_minimal() +ylab("Variant Prevalence")+
      theme(axis.text.x = element_text(angle = 45))
