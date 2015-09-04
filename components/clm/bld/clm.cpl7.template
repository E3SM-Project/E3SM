#! /bin/csh -f

if ($NINST_LND > 1) then
   set inst_counter = 1
   set inst_string = ""
   while ($inst_counter <= $NINST_LND)
      set inst_string = $inst_counter
      if ($inst_counter <= 999) set inst_string = "0$inst_string"
      if ($inst_counter <=  99) set inst_string = "0$inst_string"
      if ($inst_counter <=   9) set inst_string = "0$inst_string"
      set inst_string = _${inst_string}
      if ( ! -f "$CASEROOT/user_nl_clm${inst_string}" ) then
         cp $CODEROOT/lnd/clm/bld/user_nl_clm $CASEROOT/user_nl_clm${inst_string}
      endif
      @ inst_counter = $inst_counter + 1
   end 
else
   if ( ! -f "$CASEROOT/user_nl_clm" ) then
      cp $CODEROOT/lnd/clm/bld/user_nl_clm $CASEROOT/user_nl_clm
   endif
endif


