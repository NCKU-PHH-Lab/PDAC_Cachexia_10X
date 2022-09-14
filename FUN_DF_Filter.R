DF_Filter <- function(df,Mode = "Keyword",  # Mode = ("KeyWordFilter", "LogicalJudgment")
                      Condition = list(ColNa = "Colname", Value = 1, LogicalSet = ">"))
                     #(for Keyword) Conditio  = list(ColNa = "Colname", Exact_Match = "No", KeyWord = "Mac")
{
  if(Mode == "LogicalJudgment"){
    if(Condition[["LogicalSet"]] == ">"){
      OUTPUT <- df[df[,Condition[["ColNa"]]] > 1,]
      # OUTPUT <- df[abs(df$NES) > 1,]

    }else if(Condition[["LogicalSet"]] == ">="){
      OUTPUT <- df[df[,Condition[["ColNa"]]] >= 1,]
    }else if(Condition[["LogicalSet"]] == "<"){
      OUTPUT <- df[df[,Condition[["ColNa"]]] < 1,]
    }else if(Condition[["LogicalSet"]] == "<="){
      OUTPUT <- df[df[,Condition[["ColNa"]]] <= 1,]
    }else{
      OUTPUT <- df[df[,Condition[["ColNa"]]] == 1,]
    }

  }else{

  }





  ## Output setting ##
  OUTPUT
return(OUTPUT)


}
