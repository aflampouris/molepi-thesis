migration_events <-
  function (workdir = "~/Desktop/",
            cluster = 19,
            byGroup = "infectionRoute",
            region = "Europe",
            ntrees = 3,
            mig.type = "C",
            action = "bootstrap",
            raxml = T,
            besttree = F,
            dbuser = "R",
            schema = "thesisdb",
            dbpass = "1234") {
    set.seed(667)
    mig.type <-
      switch(mig.type,
             "C" = '==',
             "A" = "--",
             "B" = '(==|--)')
    dir.create(paste0(workdir, 'CL', cluster), showWarnings = F)
    setwd(paste0(workdir, 'CL', cluster))
    list.of.packages <-
      c(
        "parallel",
        "RMySQL",
        "tikzDevice",
        "gridExtra",
        "extrafont",
        "phangorn",
        "scales",
        "RColorBrewer",
        "phylotools",
        "reshape2",
        "ggplot2",
        "mefa",
        "plyr",
        "Hmisc",
        "stargazer",
        "papeR",
        "ips",
        "stats"
      )
    new.packages <-
      list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    if (length(new.packages))
      install.packages(new.packages)
    lapply(list.of.packages, require, character.only = T)
    filename <- paste("CL", cluster, sep = ".")
    
    #         lapply(list.of.packages,citation)
    
    
    #GET CLUSTERS#
    if (action %in% c("get.data", "bootstrap", "full", "test")) {
      dbCon <-
        dbConnect(
          MySQL(),
          user = dbuser,
          pass = dbpass,
          dbname = schema,
          host = "localhost"
        )
      dbData <- dbGetQuery(
        dbCon,
        paste0(
          "SELECT
						  uid,testYear,gender,submissionCountry,",
          byGroup,
          ",sequence
						  FROM DatasetFull d WHERE region = '",
          region,
          "' and cluster = ",
          cluster,
          ";"
        )
      )
      dbDisconnect(dbCon)
      
      dbData[, byGroup] <- factor(dbData[, byGroup], exclude = NULL)
      dbData$paupCodes <-
        factor(
          dbData[, byGroup],
          levels = unique(dbData[, byGroup]),
          exclude = NULL,
          labels = c(letters[letters != "m"], LETTERS[LETTERS != "M"], 0:9)[1:length(levels(dbData[, byGroup]))]
        )
      # dbData$paupCodes <- as.factor(dbData[["byGroup"]])
      # mapping <- data.frame('group' = unique(dbData[,byGroup]),'paupChar' = levels(dbData[["paupCodes"]]))
      mapping <- unique(dbData[, c(5, 7)])
      names(mapping) <- c("group", "paupChar")
      
      seqs <- t(sapply(strsplit(dbData[["sequence"]], ""), tolower))
      rownames(seqs) <-
        paste(dbData[["uid"]], dbData[["submissionCountry"]], dbData[, byGroup], sep = "|")
      seqs <- as.DNAbin(seqs)
      #       write.dna(seqs, file="CL23.phy", format="inter")
    }
    
    if (action %in% c("bootstrap", "full")) {
      if (raxml == T) {
        MLrun <-
          raxml(
            seqs,
            f = 'a',
            threads = 2,
            m = "GTRGAMMA",
            N = ntrees,
            p = 12345,
            x = 12345,
            k = T,
            exec = "raxmlHPC-PTHREADS-AVX"
          )
        treeObs <- mclapply(MLrun[["bootstrap"]], midpoint)
      }
      
      else if (raxml == F) {
        boot_trees <- read.tree("RAxML_bootstrap.fromR")
        treeObs <- mclapply(boot_trees, midpoint)
      }
      
      treePerm <- mclapply(treeObs, tipShuffle)
      
      write.nexus(
        treeObs,
        file = paste0(cluster, ".observed.midpoint.nex"),
        translate = F
      )
      write.nexus(
        treePerm,
        file = paste0(cluster, ".permuted.midpoint.nex"),
        translate = F
      )
      
      system(paste0("sed -i '1,/END/d' ", cluster, ".observed.midpoint.nex"))
      system(paste0("sed -i '1,/END/d' ", cluster, ".permuted.midpoint.nex"))
      
      symbols <-
        paste0("\"", paste(mapping$"paupChar", collapse = ""), "\"")
      
      PAUPhead <-
        paste0("#NEXUS \n", "BEGIN PAUP; SET MAXTREES = ", ntrees, "; END;")
      
      begin.taxa <-
        paste0(
          "BEGIN TAXA; TITLE TAXA; DIMENSIONS NTAX = ",
          dim(seqs)[1],
          "; TAXLABELS \n",
          paste(rownames(seqs), collapse = "\t"),
          "; \nEND;"
        )
      
      begin.data <-
        paste0(
          "BEGIN DATA; DIMENSIONS NTAX = ",
          dim(seqs)[1],
          " NCHAR = ",
          length(symbols),
          "; FORMAT RESPECTCASE SYMBOLS = ",
          symbols,
          " MISSING = ? GAP = - ;"
        )
      
      for (file in c("observed", "permuted")) {
        # Now actually create the PAUP input file
        nexus.filename <- paste0(filename, ".", file, ".nex")
        mesquite.filename <- paste0("mesquite.CL", cluster, ".nex")
        
        PAUPtail <- paste(
          "BEGIN PAUP;",
          "PSET OPT = DELTRAN COLLAPSE = no;",
          "SET CRITERION = parsimony TAXLABELS = full SHOWTAXNUM = yes OPT = deltran;",
          "LOG START ",
          "FILE = ",
          paste0(
            filename,
            ".",
            file,
            ".log;",
            "DESCRIBETREES 1-",
            ntrees,
            "/ BRLENS = yes CHGLIST = yes PLOT = cladogram;",
            "MPRSETS 1;LOG STOP;END;QUIT;"
          ),
          sep = '\n'
        )
        
        #MESQUITE BESTTREE
        write("#NEXUS", mesquite.filename)
        
        write(begin.data, mesquite.filename, append = T)
        write("matrix", mesquite.filename, append = T)
        write.table(
          cbind(paste0("'", rownames(seqs), "'"), as.character(dbData[["paupCodes"]])),
          mesquite.filename,
          row.names = F,
          col.names = F,
          quote = F,
          append = T
        )
        write(";\n end;", mesquite.filename, append = T)
        bestTree <- read.tree("RAxML_bestTree.fromR")
        bestTree <- midpoint(bestTree)
        
        write.nexus(bestTree, file = "RAxML_bestTree.midpoint", translate = F)
        system(paste0("sed -i '1,/END/d' RAxML_bestTree.midpoint"))
        
        file.append(mesquite.filename, "RAxML_bestTree.midpoint")
        
        #PAUP TREES
        
        write(PAUPhead, nexus.filename)
        
        if (file == "permuted") {
          write(begin.taxa, nexus.filename, append = T)
        }
        
        write(begin.data, nexus.filename, append = T)
        write("matrix", nexus.filename, append = T)
        write.table(
          cbind(paste0("'", rownames(seqs), "'"), as.character(dbData[["paupCodes"]])),
          nexus.filename,
          row.names = F,
          col.names = F,
          quote = F,
          append = T
        )
        write(";\n end;", nexus.filename, append = T)
        
        file.append(nexus.filename,
                    paste0(cluster, ".", file, ".midpoint.nex"))
        write(PAUPtail, nexus.filename, append = T)
        
        system(
          paste0(
            "wine ~/.wine/drive_c/paup95NT/win-paup4b10-console.exe -r ",
            getwd(),
            "/",
            nexus.filename
          ),
          input = c("n", "y")
        )
      }
    }
    
    # STATISTICAL TEST
    if (action %in% c("test", "full")) {
      # Define filelist
      input_type <- c("observed", "permuted")
      filelist <- paste("CL", cluster, input_type, "log", sep = ".")
      
      for (i in 1:2) {
        system(
          paste0(
            "grep -Eo '(([a-z|A-Z|0-9] ",
            mig.type,
            "> [a-z|A-Z|0-9])|Changes)' ",
            filelist[i],
            ">",
            input_type[i]
          )
        )
        
        d <-
          read.table(paste0(input_type[i]),
                     sep = " ",
                     blank.lines.skip = F)
        
        tree.start <- which(d$V1 == 'Changes')
        tree.stop <- c((tree.start - 1)[-1], nrow(d))
        if (all(tree.start - tree.stop == 0) == T) {
          next
        }
        
        get.seq <-
          mapply(
            tree.start,
            tree.stop,
            1:length(tree.start),
            FUN = function(x, y, z) {
              len <- length(seq(from = x[1], to = y[1]))
              return(rep(z, times = len))
            },
            SIMPLIFY = F
          )
        
        d$tree <- unlist(get.seq)
        d <- d[-tree.start,]
        d[, 5] <- 1
        d <- count(d, c("V1", "V2", "V3", "tree"))[,-2]
        d[, 1:2] <-
          sapply(d[, 1:2], factor, levels = mapping[["paupChar"]], labels = mapping[["group"]])
        names(d) <- c("from", "to", "bs_tree", "N")
        assign(paste0(input_type[i], ".data"), d)
        levels <- unique(d[, c(1, 2)])
        assign(paste0(input_type[i], ".levels"), levels)
      }
      
      levels <- rbind(observed.levels, permuted.levels)
      familyWise  <-
        max(length(mapping[["group"]]) * (length(mapping[["group"]]) - 1), 1)
      
      data <- merge(unique(levels), 1:ntrees)
      names(data) <- c("from", "to", "bs_tree")
      
      observed.data <-
        merge(data,
              observed.data,
              by = c("from", "to", "bs_tree"),
              all.x = T)
      observed.data$tree_type = as.factor("observed")
      permuted.data <-
        merge(data,
              permuted.data,
              by = c("from", "to", "bs_tree"),
              all.x = T)
      permuted.data$tree_type = as.factor("permuted")
      
      data <- rbind(observed.data, permuted.data)
      data[is.na(data[["N"]]), "N"] <- 0
      #     data$from = factor(data$from,levels=sort(levels(data$from),decreasing=T))
      
      data[["from"]] <-
        factor(data[["from"]],
               exclude = NULL,
               levels = sort(unique(data[["from"]]), na.last = F, decreasing = T))
      data[["to"]] <- factor(data[["to"]], exclude = NULL)
      
      paletteCols <-
        c(
          "#5E0707",
          "#6CBB1D",
          "#BB1D1D",
          "#325E07",
          "#075E5E",
          "#6C1DBB",
          "#1DBBBB",
          "#666666"
        )
      
      G3 <-
        ggplot(data, aes(x = N, fill = tree_type)) + geom_bar(alpha = 1,
                                                              position = "dodge",
                                                              width = 0.6) +
        scale_fill_manual(
          values = paletteCols[c(1, 7)],
          labels = c("Obs.", "Per."),
          name = "Tree"
        ) +
        facet_grid(from ~ to,
                   shrink = T,
                   as.table = T,
                   drop = F) +
        xlab(label = "N") + ylab(label = "Number of phylogenetic trees") +
        ggtitle("Observed Vs. permuted trees distribution") +
        theme(text = element_text(size = 20, family = "Linux Biolinum O"))
      
      
      
      data <-
        reshape(
          data,
          timevar = "tree_type",
          idvar = c("from", "to", "bs_tree"),
          direction = "wide"
        )
      #
      #       attach(data)
      #       #     data <- data[ order(data$from,data$to),]
      #       detach(data)
      #
      data <- split(data, data[, 1:2])
      data <- data[sapply(data, dim)[1,] > 0]
      
      stat_test <- as.data.frame(t(sapply(data, function(d) {
        wilcox.test(
          d$N.observed,
          d$N.permuted,
          exact = F,
          alternative = "g",
          correct = T,
          conf.int = T
        )
      })))
      stat_test <-
        cbind(colsplit(
          rownames(stat_test),
          pattern = "\\.",
          names = c("from", "to")
        ), stat_test)
      
      results <- data.frame(
        "from" = stat_test$from,
        "to" = stat_test$to,
        # 	  	"Pr." = as.numeric(stat_test$pval),
        # 	  	"Pr." = sapply(as.numeric(stat_test$pval)*familyWise ,min,0.999),
        # 	  	"est" = as.numeric(r[["estimate"]]),
        "ci" = round(as.numeric(sapply(
          stat_test[["conf.int"]], `[`, 1
        )), 1),
        #       	"Pr." = p.adjust(as.numeric(stat_test$pval),"bonferroni",n = familyWise ),
        "Pr." = p.adjust(as.numeric(stat_test$p.value), "holm", n = familyWise),
        "meanOrg" = as.numeric(lapply(sapply(data, `[`, 4), mean)),
        "meanRan" = as.numeric(lapply(sapply(data, `[`, 5), mean)),
        "medOrg" = as.numeric(lapply(sapply(data, `[`, 4), median)),
        "medRan" = as.numeric(lapply(sapply(data, `[`, 5), median))
      )
      
      
      
      
      
      results <- results[order(results[, 1]),]
      results$MeansRatio = results[["meanOrg"]] / results[["meanRan"]]
      
      
      
      N <- unique(levels)
      N1 <- as.data.frame(table(dbData[, byGroup]))
      names(N1) <- c("from", "N.from")
      N <- merge(N, N1)
      names(N1) <- c("to", "N.to")
      N <- merge(N, N1)
      N$min <- pmin(N[["N.from"]], N[["N.to"]])
      
      results <-   merge(results, N)
      results$ptr <- results[["meanOrg"]] / results$min
      results <- results[,-c(10:12)]
      
      results[["from"]] <- factor(results[["from"]], exclude = NULL)
      results[["to"]] <- factor(results[["to"]], exclude = NULL)
      
      results_tex <-
        prettify(
          results,
          extra.column = F,
          smallest.pval = -Inf,
          digits = 1,
          scientific = F
        )
      results_tex <- results_tex[,-c(1, 4)]
      results_tex$Pr. <-
        format(results[["Pr."]], scientific = T, digits = 1)
      results_tex[results_tex[["Pr."]] == '1e+00', "Pr."] <-
        "$>$0.999"
      #       results_tex[["ci"]][results$ci!=0]<-format(results_tex[["ci"]][as.numeric(results_tex$ci)>0],scientific=T,digits=0)
      # 	results_tex$ci <- format(results[["ci"]],scientific = F,digits = 3)
      # 	results_tex[["ci"]] <- paste0("[",results_tex[["ci"]],noquote(", $\\infty$)"))
      
      
      latex(
        results_tex,
        rowname = NULL,
        labels = labels(results_tex),
        size = "scriptsize",
        where = "!htpb",
        extracolsize = 'tiny'
        ,
        caption = paste(
          "Statistical testing results for monophyletic cluster No.",
          cluster
        ),
        title = paste0("tabResults", cluster)
        ,
        col.just = c(rep("l", 3), rep("c", 7)),
        na.blank = F,
        insert.bottom = paste0(attr(results_tex$"   ", "legend"), ", n=", dim(dbData)[1])
        ,
        file = paste0(
          "~/MSc-Biostatistics/Thesis/Manuscript/tables/CL",
          cluster,
          "results.tex"
        ),
        ctable = T,
        booktabs = F
        ,
        dcolumn = F,
        extracolheads = c("Group", "Group", "", "Obs.", "Perm.", "Obs.", "Perm.", "", ""),
        colheads = c(
          "Transmitter",
          "Infected",
          "$p$",
          "Mean",
          "Mean",
          "Median",
          "Median",
          "MR",
          "PTR",
          "Sig."
        )
      )
      
      write.csv(results, "results.csv")
      
      #       groupFreq <- count(dbData[,3])
      #       N <- dim(dbData)[1]
      
      results[["from"]] <-
        factor(results[["from"]],
               exclude = NULL,
               levels = sort(unique(results[["from"]]), na.last = T))
      results[["to"]] <- factor(results[["to"]], exclude = NULL)
      
      
      G1 <- ggplot(data = results, aes(x = to, y = from)) +
        geom_point(aes(color = Pr., size = MeansRatio),
                   alpha = 0.4,
                   na.rm = T) +
        geom_text(
          aes(
            label = round(ptr, 2),
            color = Pr.,
            size = MeansRatio,
            show.legend = F
          ),
          alpha = 1,
          show.legend = T
        ) +
        scale_size_continuous(
          name = "Circle=MR\nDigits=PTR",
          range = c(5, 16),
          trans = "log",
          breaks = c(0.1, 1, 2, 5, 10, 50)
        ) +
        scale_color_continuous(
          name = "p-value",
          low = paletteCols[1],
          high = paletteCols[7],
          limits = c(0, 1),
          trans = "sqrt"
        ) +
        guides(size = guide_legend()) +
        ggtitle("Tranmission events matrix") +
        theme(text = element_text(size = 20, family = "Linux Biolinum O O")) +
        ylab("Transmitting risk group") + xlab("Infected risk group")
      
      dbData$groupColors <-
        factor(
          dbData[, byGroup],
          exclude = NULL,
          levels = c("HET", "MSM", "PWID", "BTR", "OCC", "OHP", "MTC", NA),
          labels = paletteCols
        )
      groupColors <- as.character(dbData$groupColors)
      names(groupColors) <- dbData[, byGroup]
      
      dbData$genderColors <-
        factor(
          dbData[["gender"]],
          levels = c("M", "F", NA),
          exclude = NULL,
          labels = paletteCols[c(7, 6, 8)]
        )
      genderColors <- as.character(dbData$genderColors)
      names(genderColors) <- dbData[["gender"]]
      
      height <-
        18#2*max(length(levels(results[["Transmitters"]])),length(unique(country_data[["submissionCountry"]])))
      width <-
        25 # 2*max(length(levels(country_data[["submissionCountry"]])),length(unique(country_data[["submissionCountry"]])))
      
      G2 <-
        ggplot(dbData, aes_string(x = "submissionCountry", fill = byGroup)) +
        geom_bar(stat = "count", width = 0.3) + coord_flip() +
        scale_fill_manual(values = groupColors , name = "Risk\nGroup") +
        xlab(label = "Country of sampling") + ylab(label = "Number of sequences") +
        ggtitle("Submission country distribution") +
        theme(
          text = element_text(size = 20, family = "Linux Biolinum O"),
          axis.text.x = element_text(angle = 0)
        )
      
      
      G4 <-
        ggplot(dbData, aes(x = factor(testYear))) + aes_string(fill = byGroup) +
        geom_bar(stat = "count", width = 0.3) +
        scale_fill_manual(values = groupColors , name = "Risk\nGroup") +
        xlab(label = "Year of sampling") + ylab(label = "Number of sequences") +
        ggtitle("Year of sampling distribution") +
        theme(
          text = element_text(size = 20, family = "Linux Biolinum O"),
          axis.text.x = element_text(angle = 0)
        )
      
      cairo_pdf(
        "Graphs.pdf",
        width = width,
        height = height,
        family = "Linux Biolinum O",
        pointsize = 22
      )
      grid.arrange(
        G1,
        G3,
        G2,
        G4,
        ncol = 2,
        nrow = 2,
        widths = c(1, 1.2),
        heights = c(0.5, 0.5)
      )
      dev.off()
      #
    }
    return(results)
  }



for (k in c(19)) {
  migration_events(
    workdir = "~/MSc-Biostatistics/Thesis/FinalAnalysis/Testing/",
    cluster = k,
    byGroup = "iRsimple2",
    region = "Europe",
    ntrees = 50,
    mig.type = "C",
    action = "test",
    raxml = T,
    besttree = T,
    dbuser = "Ruser",
    schema = "thesisdb",
    dbpass = "1234"
  )
}