##########################################################
## This file contains Julia functions used in the Pluto ##
## notebook that replicates the results for the paper   ##
## "Identification of Incomplete Preferences" by Arie   ##
## Beresteanu and Luca Rigotti. See paper at:           ##
## https://arxiv.org/abs/2108.06282                     ##
##########################################################
## Version: 0.1                                         ##
##########################################################

function makeData(dir::String)
       # OHelections.dta was created in Stata
       file1 = string(dir,"OHelections.dta");
       dta = DataFrame(load(file1));

       #Creating a new (fixed) candidate order variable for the supreme court races.
       # The rule is that the data is sorted by county and within each county by the
       # three letter Precinct Code. Starting from the first precinct in each county
       # the candidates are ordered alphabetically. Then the order is switched (rotated)
       # after that. 
       sort!(dta, [:CountyName, :PrecinctCode]);
       dta.order =zeros(8904)

       for i=2:8904
              if dta.CountyName[i] == dta.CountyName[i-1]
                     dta.order[i]=1-dta.order[i-1]
              else
                     dta.order[i]=0
              end
       end 	


       # Read party affiliation data

		file1 = string(dir,"FullPartyAffiliation.csv");
        aff = DataFrame(CSV.File(file1));

       # The precincts are coded slightly differently in the above data so we need to fix that
       # We create a new variable in DataFrame aff which consists of two digits (01-88) and 
       # three letters. The two digits are a running count of the counties (alphabetically) and 
       # the tree letters are the three letters precinct code.
       codeTrim = s -> string(s[1:2],s[end-2:end])
       aff.code = codeTrim.(aff.PrecinctCode) 

       # We create the same variable as above for the dta DataFrame:
       sort!(dta,[:CountyName, :PrecinctCode]);
       # first creating the numerical part of the code:
       dta.codeN =ones(8904);

       for i=2:8904
              if dta.CountyName[i] == dta.CountyName[i-1]
                     dta.codeN[i]=dta.codeN[i-1]
              else
                     dta.codeN[i]=dta.codeN[i-1]+1
              end
       end
       dta.code = string.(Int.(dta.codeN),dta.PrecinctCode)
       #still not perfect becase we want county 1 to have 01 in the code (not just 1)
       shipoor = x -> length(x)==4 ? string("0",x) : x
       dta.code = shipoor.(dta.code);


       ### Joining the two data frames
       dta1 = rightjoin(aff,dta,on = :code,makeunique=true);


       ## Creating new variables for percent registered as democrat and as republican
       dta1.Dem = dta1.dem ./(dta1.dem .+ dta1.rep .+ dta1.ind)
       dta1.Rep = dta1.rep ./(dta1.dem .+ dta1.rep .+ dta1.ind)
       dta1.Ind = dta1.ind ./(dta1.dem .+ dta1.rep .+ dta1.ind)

       dta1.logTotalVoters = log.(dta1.TotalVoters);
       # adding the above didn't quite improve the regression

       ## Adding age variable
       age = DataFrame(CSV.File("c:\\data\\fullage.csv"));
       age.code = codeTrim.(age.code); #like before with aff we need to make code the same format like in dts
       dta2 = rightjoin(age,dta1,on = :code,makeunique=true);
       
       return dta2
end


function logistic(t::Real=0.0)
       return exp(t)/(1+exp(t))
end


function marginal(beta::Vector{Float64},x::Vector{Float64})
       t=exp(beta' * x);
       return beta *(t/((1+t)^2))
end


function colSum(df::DataFrame,varlist::Vector{Symbol})
	temp = combine(df, AsTable(varlist) => ByRow(sum) => :new, renamecols=false)
	return temp
end


function triLogistic(β::Vector{Float64},x1::Vector{Float64},x2::Vector{Float64})
	a = exp.([0.0, x1'*β, x2'*β]);
	return a ./ sum(a)
end

function analyzeVote(dtaCounty::DataFrame)
	#dtaCounty = dta[dta.CountyName.==c,:];
	dtaCounty.second = 1 .- dtaCounty.first;
	
	n=size(dtaCounty,1);
	

	f1 = (x,y) -> ismissing(x) ? y : x
	f2 = (x,y) -> ismissing(x) ? 1 : 0 
	
	#First Race
	varlist1 = [:CountyName, :PrecinctCode, :RegisteredVoters, :age, :Dem, :Rep, :Ind, :second, :TurnoutPercentage, :R_Baldwin];
	varlist2 = [:CountyName, :PrecinctCode, :RegisteredVoters, :age, :Dem, :Rep, :Ind, :first, :TurnoutPercentage, :R_Donnelly];
	race1 = vcat(dtaCounty[:,varlist1],dtaCounty[:,varlist2], cols=:union);
	race1.Y = f1.(race1.R_Baldwin,race1.R_Donnelly);
	race1.Donnelly = f2.(race1.R_Baldwin,race1.R_Donnelly);
	race1.First = f1.(race1.first,race1.second);
	
	#Second Race
	varlist3 = [:CountyName, :PrecinctCode, :RegisteredVoters, :age, :Dem, :Rep, :Ind,  :second, :TurnoutPercentage, :R_DeGenaro];
	varlist4 = [:CountyName, :PrecinctCode, :RegisteredVoters, :age, :Dem, :Rep, :Ind, :first, :TurnoutPercentage, :R_Stewart];
	race2 = vcat(dtaCounty[:,varlist3],dtaCounty[:,varlist4], cols=:union);
	race2.Y = f1.(race2.R_DeGenaro,race2.R_Stewart);
	race2.First = f1.(race2.first,race2.second);
	race2.Stewart = f2.(race2.R_DeGenaro,race2.R_Stewart);

	
	return n,race1,race2
end


function analyzeRaces(race1::DataFrame, race2::DataFrame)
	
	avgAge = mean(race1.ave_age);
	avgDem = mean(race1.Dem);
	avgTurnout = mean(race1.TurnoutPercentage);

	e1,e2,e3,e4,e5,e6,e7,e8 =0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
	p1,p2,p3,p4 =0.0,0.0,0.0,0.0
	
	try
		r1= reg(race1,@formula(Y ~ Donnelly+ Dem + ave_age + TurnoutPercentage +First ),Vcov.robust())
		p1=triLogistic(r1.coef,[1.0,1.0,avgDem,avgAge,avgTurnout,1.0 ],[1.0,0.0,avgDem,avgAge,avgTurnout,1.0]);
		p2=triLogistic(r1.coef,[1.0,1.0,avgDem,avgAge,avgTurnout,0.0 ],[1.0,0.0,avgDem,avgAge,avgTurnout,0.0]);
		e1,e2,e3,e4 = p1[2], p2[2],p1[3], p2[3]
	catch
		e1,e2,e3,e4 = 0,0,0,0
	end
	
	try 
		r2= reg(race2,@formula(Y ~ Stewart + Dem + ave_age + TurnoutPercentage + First ),Vcov.robust())
		p3=triLogistic(r2.coef,[1.0,1.0,avgDem,avgAge,avgTurnout,1.0 ],[1.0,0.0,avgDem,avgAge,avgTurnout,1.0]);
		p4=triLogistic(r2.coef,[1.0,1.0,avgDem,avgAge,avgTurnout,0.0 ],[1.0,0.0,avgDem,avgAge,avgTurnout,0.0]);
		e5,e6,e7,e8 = p3[2], p4[2],p3[3], p4[3]
	catch
		e5,e6,e7,e8 = 0,0,0,0
	 end
	
	p = [p1,p2,p3,p4];
	e = [e1,e2,e3,e4,e5,e6,e7,e8];
	
	return p,e
end
