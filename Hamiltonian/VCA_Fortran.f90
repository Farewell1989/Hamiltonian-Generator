subroutine gf_contract(ks,gf_buff,seqs,coords,gf_vca,ndim,ngf,ngf_vca,nclmap)
    integer,intent(in) :: ndim,ngf,ngf_vca,nclmap
    real(8),intent(in),dimension(ndim) :: ks
    complex(8),intent(in),dimension(ngf,ngf) :: gf_buff
    integer,intent(in),dimension(ngf_vca,nclmap) :: seqs
    real(8),intent(in),dimension(ngf_vca,nclmap,ndim) :: coords
    complex(8),intent(out),dimension(ngf_vca,ngf_vca) :: gf_vca
    integer :: i,j,k,l,h
    real(8),dimension(ndim) :: coords_buff
    gf_vca=(0.0_8,1.0_8)
    do i=1,ngf_vca
        do k=1,nclmap
            do j=1,ngf_vca
                do l=1,nclmap
                    do h=1,ndim
                      coords_buff(h)=coords(j,l,h)-coords(i,k,h)
                    end do
                    gf_vca(i,j)=gf_vca(i,j)+gf_buff(seqs(i,k),seqs(j,l))*exp((0.0_8,1.0_8)*dot_product(ks,coords_buff))
                end do
            end do
        end do
    end do
end subroutine gf_contract
