module analise_trelica_TCC
    
    use CargaVento
    use Estrutura_Trelica
    use LLSt
    use Matriz_Rigidez
    
    implicit none
    
    
    contains

    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S1)
    close (uni_dad_vento)
    
    
    
    
    
    
    end module