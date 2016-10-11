cluster.summary <- function( e.cutoff=1000000, nrow.cutoff=1, seq.type=names( mot.weights )[ 1 ], plot=F, sort=c("score.norm","score","resid","e.value1","e.value2","nrow") ) {
  ms <- NULL
  if ( ! is.null( seq.type ) ) ms <- meme.scores[[ seq.type ]]
  if ( is.null( ms ) ) e.cutoff <- NA
  score <-
    sapply( 1:k.clust, function( k ) mean( row.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) * resid.scaling[ iter ] +
    if ( ! is.null( mot.scores ) ) { sapply( 1:k.clust, function( k ) mean( mot.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) * mot.scaling[ iter ] } else { 0 } +
    if ( ! is.null( net.scores ) ) { sapply( 1:k.clust, function( k ) mean( net.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) * net.scaling[ iter ] } else { 0 } +
     if ( ! is.null( set.scores ) ) { sapply( 1:k.clust, function( k ) mean( set.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) * set.scaling[ iter ] } else { 0 }
  nrow <- sapply( 1:k.clust, function( k ) length( get.rows( k ) ) )

  out <- data.frame( k=1:k.clust, nrow=nrow, score=score, ##score.norm=score.norm,
                    resid=sapply( 1:k.clust, cluster.resid, varNorm=F ), 
                    consensus1=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 1 ]]$pssm ) ),
                    e.value1=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) Inf else
                      ms[[ k ]]$meme.out[[ 1 ]]$e.value ),
                    consensus2=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) "" else
                      if ( length( ms[[ k ]]$meme.out ) == 1 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 2 ]]$pssm ) ),
                    e.value2=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) Inf else
                      if ( length( ms[[ k ]]$meme.out ) <= 1 ) Inf else
                      ms[[ k ]]$meme.out[[ 2 ]]$e.value )
                    )
  if ( all( out$consensus2 == "" ) ) out$consensus2 <- out$e.value2 <- NULL
  if ( ! is.na( sort[ 1 ] ) && sort[ 1 ] %in% colnames( out ) ) out <- out[ order( out[[ sort[ 1 ] ]] ), ]
  out
}

e$cluster.summary <- cluster.summary
environment( e$cluster.summary ) <- e

