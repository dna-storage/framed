#---------------------------------------------------------------------------#
#	@Authors: Belaid Hamoum, Elsa Dupraz				    #
#	Project : DnarXiv, funded by CominLabs				    #
#	email   : belaid.hamoum@gmail.com				    #
#---------------------------------------------------------------------------#
# 			DNA Data Storage SIMULATOR			    #
#      		     (From Synthesis to Basecalling)			    #
#									    #
#---------------------------------------------------------------------------#

function channel(k,seq_set)
	nbrSim=1
	n_seq = length(seq_set)
	out_set=[]
	for strand_index=1:n_seq
	    iRef=[1]
	    startKmer=[1]
	    simSeq=[]
	    seq=seq_set[strand_index]
	    for iSeq=1:nbrSim
		nRef=length(seq)


		iRef[1]=1
		kmer=""
		#prevKmer=""
		prevEdit=""   #[1:M,2:D,3:I,4:S]

		if(k>1)
			startKmer[1]=1
		else
			startKmer[1]=2
		end



		while iRef[1] <= nRef

			currBase=seq[iRef[1]]
			if(iRef[1]>1)

				if(iRef[1]<=nRef-1)

					r1=rand()

					kmer=join(seq[startKmer[1]:iRef[1]])



					if(iRef[1]>=k )
						indexP=0

						if(get(mapKmerPrevYi[prevEdit],kmer,0)>0)
							indexP=mapKmerPrevYi[prevEdit][kmer]
						else
							indexP=mapKmerPrevYi[prevEdit]["AVG"]

						end



						#kmer INS	DEL	SUBST	ERR	MATCH
						probInsK=EditYiKmerprevYi[prevEdit][indexP,2]
						probDelK=EditYiKmerprevYi[prevEdit][indexP,3]
						probSubstK=EditYiKmerprevYi[prevEdit][indexP,4]
						probMatchK=EditYiKmerprevYi[prevEdit][indexP,6]
						probEditK=probDelK+probSubstK

						if(r1<=probEditK) #edit the kmer


								rangeDelMaxK=probDelK/probEditK
								rangeSubstMaxK=rangeDelMaxK+(probSubstK/probEditK)


								r2=rand()

								if(r2<=rangeDelMaxK) #Deletion

									deleteB(simSeq)
									prevEdit=2 #DEL
								elseif(r2<=rangeSubstMaxK) #Substitution

									substB(simSeq,rangeTransProb,currBase)
									prevEdit=4 #Subs

									r3=rand()
									if (r3<=probInsK)
											indexPLen=0
											if(get(mapKmerInsLenPrevYi[prevEdit],kmer,0)>0)
												indexPLen=mapKmerInsLenPrevYi[prevEdit][kmer]
											else
												indexPLen=mapKmerInsLenPrevYi[prevEdit]["AVG"]
											end

											rangeLenInsK=kmerInsLenProbPrevYi[prevEdit][indexPLen]
											lenInsK=length(rangeLenInsK)
											insertB(simSeq,rangeLenInsK,lenInsK)

											prevEdit=3 #INS
									end
								end


						else
							matchB(simSeq,currBase)
							prevEdit=1 #Match

							r3=rand()
							if (r3<=probInsK)
									indexPLen=0
									if(get(mapKmerInsLenPrevYi[prevEdit],kmer,0)>0)
										indexPLen=mapKmerInsLenPrevYi[prevEdit][kmer]
									else
										indexPLen=mapKmerInsLenPrevYi[prevEdit]["AVG"]
									end

									rangeLenInsK=kmerInsLenProbPrevYi[prevEdit][indexPLen]
									lenInsK=length(rangeLenInsK)
									insertB(simSeq,rangeLenInsK,lenInsK)

									prevEdit=3 #INS
							end

						end
					else

						if(seq[iRef[1]]=="A")
							if(r1<=MidProbEditA)

								r2=rand()

								if (r2<=rangeDelMaxMidA)
									deleteB(simSeq)
									prevEdit=2 #DEL
								elseif (r2<=rangeSubstMaxMidA) #substitution
									substB(simSeq,rangeTransProb,currBase)
									prevEdit=4 #SUBS

									r3=rand()
									if (r3<=MidProbInsA) #insertion
											insertB(simSeq,rangeLenInsMidA,lenInsMidA)
											prevEdit=3 #INS
									end

								end
							else
								matchB(simSeq,currBase)
								prevEdit=1 #Match

								r3=rand()
								if (r3<=MidProbInsA) #insertion
										insertB(simSeq,rangeLenInsMidA,lenInsMidA)
										prevEdit=3 #INS
								end

							end
						elseif(seq[iRef[1]]=="C")
							if(r1<=MidProbEditC)

								r2=rand()

								if (r2<=rangeDelMaxMidC)
									deleteB(simSeq)
									prevEdit=2 #DEL
								elseif (r2<=rangeSubstMaxMidC) #substitution
									substB(simSeq,rangeTransProb,currBase)
									prevEdit=4 #SUBS

									r3=rand()
									if (r3<=MidProbInsC) #insertion
											insertB(simSeq,rangeLenInsMidC,lenInsMidC)
											prevEdit=3 #INS
									end

								end
							else
								matchB(simSeq,currBase)
								prevEdit=1 #Match

								r3=rand()
								if (r3<=MidProbInsC) #insertion
										insertB(simSeq,rangeLenInsMidC,lenInsMidC)
										prevEdit=3 #INS
								end

							end
						elseif(seq[iRef[1]]=="G")
							if(r1<=MidProbEditG)

								r2=rand()


								if (r2<=rangeDelMaxMidG)
									deleteB(simSeq)
									prevEdit=2 #DEL
								elseif (r2<=rangeSubstMaxMidG) #substitution
									substB(simSeq,rangeTransProb,currBase)
									prevEdit=4 #SUBS

									r3=rand()
									if (r3<=MidProbInsG) #insertion
											insertB(simSeq,rangeLenInsMidG,lenInsMidG)
											prevEdit=3 #INS
									end

								end
							else
								matchB(simSeq,currBase)
								prevEdit=1 #Match

								r3=rand()
								if (r3<=MidProbInsG) #insertion
										insertB(simSeq,rangeLenInsMidG,lenInsMidG)
										prevEdit=3 #INS
								end

							end
						elseif(seq[iRef[1]]=="T")
							if(r1<=MidProbEditT)

								r2=rand()


								if (r2<=rangeDelMaxMidT)
									deleteB(simSeq)
									prevEdit=2 #DEL
								elseif (r2<=rangeSubstMaxMidT) #substitution
									substB(simSeq,rangeTransProb,currBase)
									prevEdit=4 #SUBS

									r3=rand()
									if (r3<=MidProbInsT) #insertion
											insertB(simSeq,rangeLenInsMidT,lenInsMidT)
											prevEdit=3 #INS
									end

								end
							else
								matchB(simSeq,currBase)
								prevEdit=1 #Match

								r3=rand()
								if (r3<=MidProbInsT) #insertion
										insertB(simSeq,rangeLenInsMidT,lenInsMidT)
										prevEdit=3 #INS
								end

							end
						end

					end
					if(iRef[1]>=k)
						#prevKmer=kmer
						startKmer[1]+=1
					end
				else



					r1=rand()

					if(r1<=EndProbEdit)

						r2=rand()

						if (r2<=rangeInsMaxEnd) #insertion
							insertB(simSeq,rangeLenInsEnd,lenInsEnd)
							prevEdit=3 #INS
						elseif (r2<=rangeDelMaxEnd) #deletion
							deleteB(simSeq)
							prevEdit=2 #DEL
						elseif (r2<=rangeSubstMaxEnd) #substitution
							substB(simSeq,rangeTransProb,currBase)
							prevEdit=4 #SUBS
						end
					else
						matchB(simSeq,currBase)
						prevEdit=1 #Match
					end
				end
			else



				r1=rand()

				if(r1<=BegProbEdit)

					r2=rand()

					if (r2<=rangeInsMaxBeg) #insertion
						insertB(simSeq,rangeLenInsBeg,lenInsBeg)
						prevEdit=3 #INS
					elseif (r2<=rangeDelMaxBeg) #deletion
						deleteB(simSeq)
						prevEdit=2 #DEL
					elseif (r2<=rangeSubstMaxBeg) #substitution
						substB(simSeq,rangeTransProb,currBase)
						prevEdit=4 #SUBS
					end
				else
					matchB(simSeq,currBase)
					prevEdit=1 #Match
				end

			end
			iRef[1]+=1;
		end
		#KV: made change so that returned values are strings, easier for pyton
		simSeq = join(simSeq)
	end
	push!(out_set,simSeq)
	end
	return out_set
end
