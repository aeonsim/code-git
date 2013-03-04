program autohmm
implicit none
integer ::i,j,k,l,io,n,npos,nclust,round,nind,ani
real*8, allocatable ::homoz(:),homozy(:,:),alpha(:,:),scaling(:),beta(:,:),gamma(:,:),ems(:,:),freq(:)
real*8, allocatable ::outprob(:,:),num_bjk(:,:),den_bjk(:),num_pi(:),den_pi(:),trans(:,:)
real*8 ::f0,f1,f2,ptrans,pinit(3),lik,val,num_transition

nind=25

open(10,file='homozygosity_id.txt') ! file with homozygosity at each marker / or probability of homozygosity
open(11,file='mafs.txt')

npos=0
do
read(11,*,iostat=io)i,f0
if(io/=0)exit
npos=npos+1
enddo

allocate(homozy(npos,nind),homoz(npos),freq(npos),outprob(npos,nind))

rewind(11)

n=0;freq=0.00
do
 read(11,*,iostat=io)i,f1
 if(io/=0)exit
 n=n+1
 freq(n)=f1
enddo

n=1;homozy=0.00
do
 read(10,*,iostat=io)i,(homozy(n,j),j=1,25)
 if(io/=0)exit
 n=n+1
! homoz(n)=f1
enddo
print*,'Number of positions red ::',n

!############ INITIALIZE PARAMETERS ################

! K, number of states

nclust=2

! transition probabilities

allocate(trans(nclust,nclust))
ptrans=0.0001 ! f(H1)/f(H0) approximately ptrans**2 to go into state
trans=ptrans
trans(1,1)=1.0-ptrans
trans(2,2)=1.0-ptrans

pinit(1)=0.99999;pinit(2)=0.00001

! emission probabilities => see function

!###### start HMM

allocate(alpha(nclust,npos),scaling(npos))
allocate(beta(nclust,npos),gamma(nclust,npos))

do ani=1,25
 homoz=homozy(:,ani)

!############ FORWARD ALGORITHM ####################

alpha=0.d0;scaling=0.d0

! initialisation

do i=1,nclust
  alpha(i,1)=pinit(i)*emission(i,1)
  scaling(1)=scaling(1)+alpha(i,1)
enddo

 scaling(1)=1.0/scaling(1)
 alpha(:,1)=alpha(:,1)*scaling(1)

! induction: to get to i, two ways:
! 1) transition + jump into cluster i
! 2) no transition and in i at previous position
 
 do k=2,npos
  do i=1,nclust
   do j=1,nclust
    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*trans(j,i)*emission(i,k)
   enddo
    scaling(k)=scaling(k)+alpha(i,k)
  enddo
 scaling(k)=1.0/scaling(k)
 alpha(:,k)=alpha(:,k)*scaling(k)
 enddo  

! termination

 lik=-sum(log(scaling(:)))
 print*,ani,lik

!############ BACKWARD ALGORITHM ####################

beta=0.d0;gamma=0.d0

! initialisation

do i=1,nclust
 gamma(i,npos)=alpha(i,npos)*1.0  ! beta(i,j,npos)=1.0
 beta(i,npos)=1.0*scaling(npos)
enddo

! induction
! to arrive in k: with or without transition

do k=npos-1,1,-1
 do i=1,nclust
  do j=1,nclust
    beta(i,k)=beta(i,k)+trans(i,j)*emission(j,k+1)*beta(j,k+1)
  enddo
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)
 enddo
enddo

do k=1,npos
 gamma(:,k)=gamma(:,k)/sum(gamma(:,k))
enddo

outprob(1:npos,ani)=gamma(2,:)

enddo ! ani


! ****** OUTPUT *********

open(12,file='autozygosity.prob')
do i=1,npos
! write(12,'(i7,<nind>(1x,f20.17))')i,gamma(1,i),gamma(2,i)
 write(12,'(i7,25(1x,f20.17))')i,(outprob(i,j),j=1,nind)
enddo

contains

!########### EMISSION: NON-PARAMETRIC ################

function emission(cluster1,marker)
implicit none
integer ::cluster1,marker
real*8 ::emission,maf,f1,f2

f1=freq(marker);f2=1.0-freq(marker)
maf=min(f1,f2)

if(maf>0.10 .and. homoz(marker)/=2.00)then ! use marker
 if(cluster1==1)then
   emission=homoz(marker)*(1.0-2.0*f1*f2)+(1.00-homoz(marker))*(2*f1*f2)
 else if(cluster1==2)then
   emission=homoz(marker)*1.000+(1.00-homoz(marker))*0.001
 endif
else ! skip marker with low quality or maf
 emission=1.00
endif

end function

end program
