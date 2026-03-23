########### This Module was written by Dr. Jules Marien, Postdoctoral researcher, Laboratoire de Biochimie Theorique (2026) ################

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt





def contact_map_atom2residue(contact_map, groupe_1_indices_change, groupe_2_indices_change, groupe_1_resnames, groupe_2_resnames):
    """MDAnalysis can only compute contact between atoms, but does not define the contact between residues. This function transforms an atomic contact map into a residue contact map by using the indices of change of residues. A contact between residues is therefore defined as any atom in residue 1 being at less than cutoff to any atom in residue 2. We advise to remove hydrogens from selections in Contact_map_calculation to get more interesting and faster calculations. """

    #Initialize the residue contact matrix as a boolean matrix filled with False
    contact_matrix_residues = np.zeros((len(groupe_1_resnames),  len(groupe_2_resnames) ), dtype=bool)

    #Loop over residues of groupe_1 
    for i in range(1, len(groupe_1_resnames)+1):
        index_atom_start_residue_groupe_1 = groupe_1_indices_change[i-1] 
        index_atom_last_residue_groupe_1 = groupe_1_indices_change[i] 


        #Loop over residues of groupe_2 
        for j in range(1, len(groupe_2_resnames)+1):
            index_atom_start_residue_groupe_2 = groupe_2_indices_change[j-1] 
            index_atom_last_residue_groupe_2 = groupe_2_indices_change[j] 
            
            #The next line tells if there's a contact happening in the submatrix defined as (atoms of residue of groupe_1) x (atoms of residue of groupe_2)
            contact_per_residue = np.any( contact_map[index_atom_start_residue_groupe_1:index_atom_last_residue_groupe_1+1, index_atom_start_residue_groupe_2:index_atom_last_residue_groupe_2+1]  )

            #If there is a contact between residues, we put True
            if contact_per_residue == True:
                contact_matrix_residues[i-1,j-1] = True

    return contact_matrix_residues






def Contact_map_calculation(coordinates_file_name, trajectory_file_name, selection_1, selection_2, cutoff, name_and_location_contact_map_atomic, name_and_location_contact_map_residue, name_and_location_number_of_contacts_atomic, name_and_location_number_of_contacts_residue, name_and_location_residues_number_and_name_groupe_1, name_and_location_residues_number_and_name_groupe_2, name_and_location_atom_number_and_name_groupe_1, name_and_location_atom_number_and_name_groupe_2):


    #Check if the selections are the same. We will need to remove the upper matrix and the diagonal if it's the case. 
    #/!\ The numbers of contacts must not be trusted if your selections overlap ! We only took care of the case when the selections are the same
    same_selection = False
    if selection_1 == selection_2:
        same_selection = True

    #Load the trajectory
    u = mda.Universe(coordinates_file_name, trajectory_file_name)

    #Make the atom selections
    u.trajectory[:]
    groupe_1 = u.select_atoms( selection_1)

    u.trajectory[:]
    groupe_2 = u.select_atoms( selection_2)

    #Get the resid numbers and residue names per atom

    groupe_1_resid_numbers_per_atom = groupe_1.resnums
    groupe_1_resnames_per_atom = groupe_1.resnames


    groupe_2_resid_numbers_per_atom = groupe_2.resnums
    groupe_2_resnames_per_atom = groupe_2.resnames


    #Create a list containing the indices for which the residue change

    groupe_1_indices_change = [0]
    for i in range(1,len(groupe_1_resnames_per_atom)):
        if groupe_1_resid_numbers_per_atom[i-1] != groupe_1_resid_numbers_per_atom[i]:
            groupe_1_indices_change.append(i)
    #Add the last index by hand
    groupe_1_indices_change.append(len(groupe_1_resid_numbers_per_atom))
   


    groupe_2_indices_change = [0]
    for i in range(1,len(groupe_2_resnames_per_atom)):
        if groupe_2_resid_numbers_per_atom[i-1] != groupe_2_resid_numbers_per_atom[i]:
            groupe_2_indices_change.append(i)
    #Add the last index by hand
    groupe_2_indices_change.append(len(groupe_2_resid_numbers_per_atom))


    #Get the residue names list (per residue)
    groupe_1_resnames_per_residue = groupe_1.residues.resnames
    groupe_2_resnames_per_residue = groupe_2.residues.resnames

    #We can use the indices in groupe_1_indices_change and groupe_1_indices_change to target the proper elements in the atomic contact matrix !

    
    #initialize the calculation by calculating the contact map of frame 0 and 1
    u.trajectory[0]
    distance_matrix_per_atom = contacts.distance_array(groupe_1.positions, groupe_2.positions)

    # determine which distances <= radius
    cumulative_matrix_atomic_contacts = contacts.contact_matrix(distance_matrix_per_atom, cutoff)

    #Transform into residue contacts
    cumulative_matrix_residue_contacts = contact_map_atom2residue(cumulative_matrix_atomic_contacts, groupe_1_indices_change, groupe_2_indices_change, groupe_1_resnames_per_residue, groupe_2_resnames_per_residue)


    #Check if the selections overlap
    if same_selection == False:

        #Number of atomic contacts over time
        number_atomic_contacts_over_time = [cumulative_matrix_atomic_contacts.astype(int).sum()]

        #Number of residue contacts over time
        number_residue_contacts_over_time = [cumulative_matrix_residue_contacts.astype(int).sum()]

    #If the selections are the same, we must remove the upper matrix and the diagonal (the diagonal will be fully positive since by definition an atom is in contact with itself all the time)
    else:

        #Number of atomic contacts over time
        number_atomic_contacts_over_time = [np.tril(cumulative_matrix_atomic_contacts.astype(int)).sum() - len(cumulative_matrix_atomic_contacts[0,:])]

        #Number of residue contacts over time
        number_residue_contacts_over_time = [np.tril(cumulative_matrix_residue_contacts.astype(int)).sum() - len(cumulative_matrix_residue_contacts[0,:])]  
    


    #We calculate for each frame 
    #To gain memory, we add the calculation at each step so that we only have 1 matrix of size groupe_1 x groupe_2. We just need to divide by the number of frames and multiply by 100 at the end to get the percentage of contact time

    #We convert the boolean arrays into integer arrays

    cumulative_matrix_atomic_contacts.astype(int)
    cumulative_matrix_residue_contacts.astype(int)

    u.trajectory[:]

    number_of_frames = len(u.trajectory)

    for ts in u.trajectory[1:]:
        #groupe_1.positions and groupe_2.positions change as ts changes  
        distance_matrix_per_atom = contacts.distance_array(groupe_1.positions, groupe_2.positions) 

        matrix_atomic_contacts_at_ts = contacts.contact_matrix(distance_matrix_per_atom, cutoff)

        matrix_residue_contacts_at_ts = contact_map_atom2residue(matrix_atomic_contacts_at_ts, groupe_1_indices_change, groupe_2_indices_change, groupe_1_resnames_per_residue, groupe_2_resnames_per_residue)

        #print(matrix_residue_contacts_at_ts)

        cumulative_matrix_atomic_contacts = cumulative_matrix_atomic_contacts + matrix_atomic_contacts_at_ts.astype(int)

        cumulative_matrix_residue_contacts = cumulative_matrix_residue_contacts + matrix_residue_contacts_at_ts.astype(int)

        #Check if the selections overlap
        if same_selection == False:
            number_atomic_contacts_over_time.append(matrix_atomic_contacts_at_ts.astype(int).sum())

            number_residue_contacts_over_time.append(matrix_residue_contacts_at_ts.astype(int).sum())

        #If the selections are the same, we must remove the upper matrix and the diagonal (the diagonal will be fully positive since by definition an atom is in contact with itself all the time)
        else:
            number_atomic_contacts_over_time.append(np.tril(matrix_atomic_contacts_at_ts.astype(int)).sum() - len(matrix_atomic_contacts_at_ts[0,:]))

            number_residue_contacts_over_time.append(np.tril(matrix_residue_contacts_at_ts.astype(int)).sum() - len(matrix_residue_contacts_at_ts[0,:]))

    #Divide the matrices by the number of frames  and multiply by 100 to get the percentage of contact 

    matrix_atomic_contacts = cumulative_matrix_atomic_contacts / number_of_frames *100

    matrix_residue_contacts = cumulative_matrix_residue_contacts /number_of_frames *100


    #Get the resid numbers uniquely 

    groupe_1_resid_numbers_per_residue = groupe_1.residues.resnums

    groupe_2_resid_numbers_per_residue = groupe_2.residues.resnums



    #Save the results
    np.savetxt(name_and_location_contact_map_atomic+".txt", matrix_atomic_contacts, header="Shape of the matrix : "+str(len(groupe_1))+" x "+str(len(groupe_2)))

    np.savetxt(name_and_location_contact_map_residue+".txt", matrix_residue_contacts, header="Shape of the matrix : "+str(len(groupe_1_resnames_per_residue))+" x "+str(len(groupe_2_resnames_per_residue)))

    with open(name_and_location_number_of_contacts_atomic+".txt", "w") as f:
        for items in number_atomic_contacts_over_time:
            f.write('%i\n' %items)

    with open(name_and_location_number_of_contacts_residue+".txt", "w") as f:
        for items in number_residue_contacts_over_time:
            f.write('%i\n' %items)


    #Save the residue numbers and names for each group for the residue contact map

    residues_number_and_name_groupe_1 = []
    for i in range(len(groupe_1_resnames_per_residue)):
        residues_number_and_name_groupe_1.append(str(groupe_1_resid_numbers_per_residue[i]) + " " + groupe_1_resnames_per_residue[i])

    residues_number_and_name_groupe_2 = []
    for i in range(len(groupe_2_resnames_per_residue)):
        residues_number_and_name_groupe_2.append(str(groupe_2_resid_numbers_per_residue[i]) + " " + groupe_2_resnames_per_residue[i])

    np.savetxt(name_and_location_residues_number_and_name_groupe_1+".txt", residues_number_and_name_groupe_1, fmt="%s")

    np.savetxt(name_and_location_residues_number_and_name_groupe_2+".txt", residues_number_and_name_groupe_2, fmt="%s")


    #Save the residue numbers and names and atom type for each group for the atomic contact map

    atom_number_and_name_groupe_1 = []
    for i in range(len(groupe_1.resnames)):
        atom_number_and_name_groupe_1.append(str(groupe_1.atoms.resnums[i]) + " " + str(groupe_1.atoms.resnames[i]) + " "+ str(groupe_1.atoms.names[i]))

    atom_number_and_name_groupe_2 = []
    for i in range(len(groupe_2.resnames)):
        atom_number_and_name_groupe_2.append(str(groupe_2.atoms.resnums[i]) + " " + str(groupe_2.atoms.resnames[i]) + " "+ str(groupe_2.atoms.names[i]))


    np.savetxt(name_and_location_atom_number_and_name_groupe_1+".txt", atom_number_and_name_groupe_1, fmt="%s")

    np.savetxt(name_and_location_atom_number_and_name_groupe_2+".txt", atom_number_and_name_groupe_2, fmt="%s")








def Filter_contact_map(file_contact_map, file_names_groupe_1, file_names_groupe_2, percentage_cutoff, name_and_location_filtered_contact_map, name_and_location_filtered_list_names_groupe_1, name_and_location_filtered_list_names_groupe_2):

    #Get the shape of the matrix thanks to the header
    with open(file_contact_map, "r") as f:
        header = f.readline()

    dimension_x = int( header.split()[-3] )
    dimension_y = int( header.split()[-1] )
    print(dimension_x)
    print(dimension_y)

    #Load the data
    contact_map = np.loadtxt(file_contact_map)
    contact_map.reshape(dimension_x, dimension_y)

    list_names_file_groupe1 = open(file_names_groupe_1)
    list_names_groupe_1 = list_names_file_groupe1.readlines()

    list_names_file_groupe2 = open(file_names_groupe_2)
    list_names_groupe_2 = list_names_file_groupe2.readlines()


    #Create an array with True where there's contact above percentage_cutoff and False otherwise, with the same size as contact_map

    filter_array_along_y = np.empty(np.shape(contact_map),dtype=bool)

    for i in range(len(contact_map[:,0])):
    
        for j in range(len(contact_map[0,:])):
        
            if contact_map[i,j] > percentage_cutoff:
            
                filter_array_along_y[i,j] = True
            
            else:
                filter_array_along_y[i,j] = False



    #return the list of indices in the contact map of the residues in groupe_2 being in contact with groupe_1

    len_groupe_1 = len(list_names_groupe_1)
    len_groupe_2 = len(list_names_groupe_2)

    residue_groupe_1 = 0

    list_indices_contacts_along_y = []

    for j in range(len_groupe_2):
    
        while filter_array_along_y[residue_groupe_1,j] == False:
        
            residue_groupe_1 = residue_groupe_1 +1
        
            if residue_groupe_1 == len_groupe_1 :
                break
        
        if residue_groupe_1 < len_groupe_1:
        
            list_indices_contacts_along_y.append(j)
        
        residue_groupe_1 = 0



    #Recreating the filtered matrix in the y dimension             
                
    contact_map_filtered_along_y = np.empty((len_groupe_1,len(list_indices_contacts_along_y)))

    for j in range(len(list_indices_contacts_along_y)):

        contact_map_filtered_along_y[:,j] = contact_map[:,list_indices_contacts_along_y[j]] 


    #Filtering the list of names along y
    filtered_names_contacts_along_y = []

    for j in range(len(list_indices_contacts_along_y)):
    
        filtered_names_contacts_along_y.append(list_names_groupe_2[list_indices_contacts_along_y[j]])


    #Repeat along x :

    #Create an array with True where there's contact above percentage_cutoff and False otherwise, with the same size as contact_map

    filter_array_along_x = np.empty(np.shape(contact_map_filtered_along_y),dtype=bool)

    for i in range(len(contact_map_filtered_along_y[:,0])):
    
        for j in range(len(contact_map_filtered_along_y[0,:])):
        
            if contact_map_filtered_along_y[i,j] > percentage_cutoff:
            
                filter_array_along_x[i,j] = True
            
            else:
                filter_array_along_x[i,j] = False



    #return the list of indices in the contact map of the residues in groupe_2 being in contact with groupe_1

    residue_groupe_2 = 0

    list_indices_contacts_along_x = []

    for i in range(len_groupe_1):
    
        while filter_array_along_x[i,residue_groupe_2] == False:
        
            residue_groupe_2 = residue_groupe_2 +1
        
            if residue_groupe_2 == len(list_indices_contacts_along_y) :
                break
        
        if residue_groupe_2 < len(list_indices_contacts_along_y):
        
            list_indices_contacts_along_x.append(i)
        
        residue_groupe_2 = 0



    #Recreating the filtered matrix in both dimensions             
                
    contact_map_filtered = np.empty((len(list_indices_contacts_along_x),len(list_indices_contacts_along_y)))

    for i in range(len(list_indices_contacts_along_x)):

        contact_map_filtered[i,:] = contact_map_filtered_along_y[list_indices_contacts_along_x[i],:] 


    #Filtering the list of names along x
    filtered_names_contacts_along_x = []

    for i in range(len(list_indices_contacts_along_x)):
    
        filtered_names_contacts_along_x.append(list_names_groupe_1[list_indices_contacts_along_x[i]])


    #Remove the "\n" character at the end of each string in the lists

    treated_filtered_names_contacts_along_x = []
    for elem in filtered_names_contacts_along_x:
        treated_filtered_names_contacts_along_x.append(elem[:-1])

    treated_filtered_names_contacts_along_y = []
    for elem in filtered_names_contacts_along_y:
        treated_filtered_names_contacts_along_y.append(elem[:-1])


    #Save the data

    np.savetxt(name_and_location_filtered_contact_map+".txt", contact_map_filtered, header="Shape of the matrix : "+str(len(filtered_names_contacts_along_x))+" x "+str(len(filtered_names_contacts_along_y)))

    np.savetxt(name_and_location_filtered_list_names_groupe_1+".txt", treated_filtered_names_contacts_along_x, fmt="%s")

    np.savetxt(name_and_location_filtered_list_names_groupe_2+".txt", treated_filtered_names_contacts_along_y, fmt="%s")














import matplotlib.ticker as mticker

def Plot_contact_map(file_contact_map, file_names_groupe_1, file_names_groupe_2, size_figure_tuple, grid_boolean, cmap, v_min, v_max, fontsize_colorbar, xlabel, ylabel, rotation, name_and_location, colorbar_tick_fontsize, axis_tick_fontsize):

    #Get the shape of the matrix thanks to the header
    with open(file_contact_map, "r") as f:
        header = f.readline()

    dimension_x_contact_map = int( header.split()[-3] )
    dimension_y_contact_map = int( header.split()[-1] )
    print(dimension_x_contact_map)
    print(dimension_y_contact_map)

    #Load the data
    contact_map = np.loadtxt(file_contact_map)
    contact_map.reshape(dimension_x_contact_map, dimension_y_contact_map)

    list_names_file_groupe1 = open(file_names_groupe_1)
    #list_names_groupe_1 = list_names_file_groupe1.readlines()

    list_names_groupe_1 = [elem.strip('\n') for elem in list_names_file_groupe1.readlines()]

    list_names_file_groupe2 = open(file_names_groupe_2)
    #list_names_groupe_2 = list_names_file_groupe2.readlines()

    list_names_groupe_2 = [elem.strip('\n') for elem in list_names_file_groupe2.readlines()]

    #Plot
    plt.figure(figsize=size_figure_tuple, dpi=600)
    ax=plt.axes()

    plt.grid(grid_boolean)

    #plt.xlim(-0.5,len(list_names_groupe_1)-0.5)
    #plt.ylim(-0.5,len(list_names_groupe_2)-0.5)

    #Careful to the dimensions here. I don't really know why we need to transpose the matrix but that's how it works 
    #plt.imshow(contact_map.T, cmap=cmap,vmin =v_min, vmax =v_max)
    plt.imshow(contact_map, cmap=cmap,vmin =v_min, vmax =v_max)

    cbar = plt.colorbar()
    cbar.set_label(label='% of contacts', fontsize=fontsize_colorbar, weight='bold')
    cbar.ax.tick_params(labelsize=colorbar_tick_fontsize)

    #Setting ticks 

    if dimension_x_contact_map < dimension_y_contact_map: 
        ax.xaxis.set_major_locator(mticker.FixedLocator(range(len(list_names_groupe_1))))
        ax.xaxis.set_major_formatter(mticker.FixedFormatter(list_names_groupe_1))

        ax.yaxis.set_major_locator(mticker.FixedLocator(range(len(list_names_groupe_2))))
        ax.yaxis.set_major_formatter(mticker.FixedFormatter(list_names_groupe_2))

        plt.xlim(-0.5,len(list_names_groupe_1)-0.5)
        plt.ylim(-0.5,len(list_names_groupe_2)-0.5)

    else :
        ax.xaxis.set_major_locator(mticker.FixedLocator(range(len(list_names_groupe_2))))
        ax.xaxis.set_major_formatter(mticker.FixedFormatter(list_names_groupe_2))

        ax.yaxis.set_major_locator(mticker.FixedLocator(range(len(list_names_groupe_1))))
        ax.yaxis.set_major_formatter(mticker.FixedFormatter(list_names_groupe_1))

        plt.xlim(-0.5,len(list_names_groupe_2)-0.5)
        plt.ylim(-0.5,len(list_names_groupe_1)-0.5)
    

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.xticks(rotation=rotation, fontsize=axis_tick_fontsize)
    plt.yticks(fontsize=axis_tick_fontsize)

    plt.tight_layout()

    plt.savefig(name_and_location+'.png')



