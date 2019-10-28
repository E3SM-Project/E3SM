      subroutine five_wrapper(pver_five_out)
          use five_intr, only: pver_five

          integer pver_five_out 

          pver_five_out = pver_five

      end subroutine five_wrapper 
