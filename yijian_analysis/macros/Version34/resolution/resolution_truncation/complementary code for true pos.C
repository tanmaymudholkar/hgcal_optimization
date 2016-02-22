      std::vector<HGCSSGenParticle> * genvec = 0; // get the generated particle vector 

      lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec); 

        if((*genvec).size() > 1) continue;  // Delete events which has converted photons.
         else{   
             const HGCSSGenParticle gHit = (*genvec)[0];
    
	     double gposx = gHit.x();
	     double gposy = gHit.y();
             double geta = gHit.eta();
             double gphi = gHit.phi();
             }

