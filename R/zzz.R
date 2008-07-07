#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

.First.lib <- function(lib, pkg) {
    library.dynam("caTools", pkg, lib)
}
