	variable tcycleLG 0.816497e6/$2
	variable PMLGdt $tcycleLG/$1
	variable PMLGinc 207.8/$1
	pulse $PMLGdt $2 [0.5*$PMLGinc:$PMLGinc:207.8-0.5*$PMLGinc]
	pulse $PMLGdt $2 [27.8-0.5*$PMLGinc:-$PMLGinc:-180+0.5*$PMLGinc]
