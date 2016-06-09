Imports System.IO

Public Class clsIOData
    Public Const STR_SIZE As Short = 14
    Public Const COMMENT_SIZE As Short = 80
    Public Const EPC9_STATS_SIZE As Short = 104
    Public Const ROOT_TEXT_SIZE As Short = 400
    Public Const MACRO_NAME_SIZE As Short = 16

    ' rr   - root record
    ' tr   - tree record

    ' srr  - series record
    ' swr  - sweep_record

    ' smr  - stimulation_record
    ' sgr  - segment_record

    Structure TreeRecord
        Dim magic_number As Int32
        Dim number_tree_levels As Int32
        Dim tree_level() As Int32
        Public Sub RDim(ByVal i As Int32)
            ReDim tree_level(i)
        End Sub
    End Structure

    ' .pul structures
    Structure PulRootRecord
        Dim version As Int32
        Dim version_name14 As String
        Dim file_name14 As String
        Dim comments400 As String
        Dim start_time As Double
        Dim number_children As Int32
        Dim gr() As PulGroupRecord
        Public Sub Rdim(ByVal i As Int32)
            ReDim gr(i)
        End Sub
    End Structure

    Structure PulGroupRecord
        Dim label14 As String
        Dim text80 As String
        Dim experimenter_number As Int32
        Dim extra_long_real As Double

        Dim number_children As Int32
        Dim srr() As PulSeriesRecord
        Public Sub Rdim(ByVal i As Int32)
            ReDim srr(i)
        End Sub
    End Structure

    Structure PulSeriesRecord
        Dim time As Double
        Dim bandwidth As Double
        Dim pipette_potential As Double
        Dim cell_potential As Double
        Dim pipette_resistance As Double
        Dim seal_resistance As Double
        Dim background_noise As Double
        Dim temperature As Double
        Dim pipette_pressure As Double
        Dim user_param1_value As Double
        Dim user_param1_name14 As String
        Dim user_param1_unit As Int16
        Dim user_param2_value As Double
        Dim user_param2_name14 As String
        Dim user_param2_unit As Int16
        Dim recording_mode As Byte
        Dim filler1 As Byte
        Dim comment80 As String
        Dim epc9_state104 As String
        Dim internal_solution As Int32
        Dim external_solution As Int32
        Dim extra_y_unit1 As Int16
        Dim extra_y_unit2 As Int16
        Dim disp_y_unit1 As Int32
        Dim disp_y_unit2 As Int32
        Dim fura_k As Double
        Dim fura_min As Double
        Dim fura_max As Double
        Dim lock_in_ext_phase As Double
        Dim timer As Double
        Dim extra_long_real As Double
        '................
        Dim swr() As PulSweepRecord
        Dim number_children As Int32
        Public Sub Rdim(ByVal i As Int32)
            ReDim swr(i)
        End Sub
    End Structure

    Structure PulSweepRecord
        Dim time As Double
        Dim stim_count As Int32
        Dim sweep_count As Int32
        Dim average_count As Int32
        Dim leak As Byte
        Dim second_trace As Byte
        Dim label14 As String
        Dim data_points As Int32
        Dim data As Int32
        Dim data_pointer As Int32
        Dim data_factor1 As Double
        Dim data_factor2 As Double
        Dim c_slow As Double
        Dim g_series As Double
        Dim rs_value As Double
        Dim rs_fraction As Double
        Dim zero_current As Double
        Dim online_y_result As Double
        Dim online_x_result As Double
        Dim total_points As Int32
        Dim offset As Int32
        Dim sweep_kind As Int16
        Dim fura_points As Int32
        Dim fura_data As Int32
        Dim fura_pointer As Int32
        Dim online_y_result2 As Double
        Dim online_x_result2 As Double
        Dim disp_factor1 As Double
        Dim disp_factor2 As Double
        Dim data_format As Byte
        Dim data_abscissa As Byte
        Dim timer As Double
        Dim sweep_spares As Int32
        '........................
        Dim number_children As Int32   '  = 0  now
    End Structure
    ' .pgf file

    Structure PgfRootRecord
        Dim root_version As Int16
        Dim number_children As Int32
        Dim smr() As PgfStimulationRecord
        Public Sub RDim(ByVal i As Int32)
            ReDim smr(i)
        End Sub
    End Structure

    Structure PgfStimulationRecord
        Dim file_name14 As String
        Dim entry_name14 As String
        Dim sample_interval As Double     ' sampling rate period
        Dim filter_factor As Double
        Dim st_sweep_interval As Double
        Dim number_sweeps As Int32
        Dim number_repeats As Int32
        Dim repeat_wait As Double
        Dim linked_sequence14 As String
        Dim linked_wait As Double
        Dim leak_count As Int32
        Dim leak_size As Double
        Dim leak_holding As Double
        Dim leak_alternate As Byte
        Dim alt_leak_averaging As Byte
        Dim leak_delay As Double
        Dim trigger_segment1 As Int16
        Dim trigger_time1 As Double
        Dim trigger_length1 As Double
        Dim trigger_amplitude1 As Double
        Dim trigger_dac1 As Int16
        Dim trigger_segment2 As Int16
        Dim trigger_time2 As Double
        Dim trigger_length2 As Double
        Dim trigger_amplitude2 As Double
        Dim trigger_dac2 As Int16
        Dim trigger_segment3 As Int16
        Dim trigger_time3 As Double
        Dim trigger_length3 As Double
        Dim trigger_amplitude3 As Double
        Dim trigger_dac3 As Int16
        Dim number_of_triggers As Int16
        Dim relevant_x_segment As Int16
        Dim relevant_y_segment As Int16
        Dim write_mode As Byte
        Dim increment_mode As Byte
        Dim total_sweep_length As Int32
        Dim max_sweep_length As Int32
        Dim input_channels As Int16
        Dim g_update As Byte
        Dim rel_abs_pot As Byte
        Dim has_continuous As Byte
        Dim log_increment As Byte
        Dim pgf_stim_dac As Int16
        Dim adc1 As Int16
        Dim adc2 As Int16
        Dim y_unit1 As Int16
        Dim y_unit2 As Int16
        Dim v_memb_increment As Int32
        Dim ext_trigger As Byte
        Dim file_template As Byte
        Dim stim_kind As Int16
        Dim lock_in_cycle As Double
        Dim lock_in_amplitude As Double
        Dim fura_on As Byte
        Dim v_memb_mode As Byte     'pgf_st_vmemb_mode
        Dim fura_tot_length As Double
        Dim fura_delay As Double
        Dim fura_length1 As Double
        Dim fura_length2 As Double
        Dim fura_wave_length0 As Double
        Dim fura_wave_length1 As Double
        Dim fura_wave_length2 As Double
        Dim fura_repeats As Int16
        Dim lock_in_skip As Int32
        Dim lock_in_v_reversal As Double
        Dim lock_in_mode As Byte
        Dim lock_in_show As Byte
        Dim config_macro16 As String
        Dim end_macro16 As String
        Dim ampl_mode_kind As Byte
        Dim filler1 As Byte
        '...........
        Dim number_children As Int32
        Dim sgr() As PgfSegmentRecord
        Public Sub Rdim(ByVal i As Int32)
            ReDim sgr(i)
        End Sub
    End Structure

    Structure PgfSegmentRecord
        Dim se_class As Byte
        Dim is_holding As Byte
        Dim voltage As Double
        Dim duration As Double
        Dim delta_v_factor As Double
        Dim delta_v_increment As Double
        Dim delta_t_factor As Double
        Dim delta_t_increment As Double
        ' ..........
        Dim number_children As Int32
    End Structure

    Structure IData
        Dim time As Double
        Dim trace1 As Double  ' signal
        Dim trace2 As Double  ' data
    End Structure

    Structure InData
        Dim sweep As String
        Dim rate As Double
        Dim size As Integer ' the actual size of id array
        Dim fname As String
        Dim id() As IData
        Public Sub New(ByVal index As Integer)
            ReDim id(index)
            size = 0
        End Sub
        Public Overloads Sub Copy(ByVal time As Double, ByVal trace1 As Double, ByVal trace2 As Double, ByVal i As Integer)
            id(i).time = time       ' in ms
            id(i).trace1 = trace1
            id(i).trace2 = trace2
        End Sub

        Public Overloads Sub Copy(ByVal time As Double, ByVal trace2 As Double, ByVal i As Integer)
            id(i).time = time
            id(i).trace2 = trace2
        End Sub
    End Structure
    Public org_xmin_value As Double
    Public org_xmax_value As Double
    Public cur_xmin_value As Double
    Public cur_xmax_value As Double

    Public org_ymax_value As Double
    Public org_ymin_value As Double
    Public cur_ymax_value As Double
    Public cur_ymin_value As Double

    Private buffer(1000) As Byte
    Private pul_tr, pgf_tr As TreeRecord
    Private pul_rr As PulRootRecord
    Private pgf_rr As PgfRootRecord

    Public in_data(MAXSWEEPINDEX) As InData
    Private p As ProgressBar
    Private bytes As Integer
    'Private t As Thread
    Public sizei As Short


    Public Function ReadDatFile(ByVal dir_name As String, ByVal n As Short) As Short

        Dim i, j, k, l As Integer
        Dim fname As String = Nothing

        Dim ns As Short = 0
        Dim tmp As Double
        ' Exit if the inpit file hasn't been choosen
        If (dir_name = "") Then
            Return ns
            Exit Function
        End If

        org_ymax_value = -9999
        org_ymin_value = 9999
        cur_ymax_value = -9999
        cur_ymin_value = 9999


        ' get file name without extention
        fname = dir_name.Remove(dir_name.LastIndexOf("."), dir_name.Length - dir_name.LastIndexOf("."))

        If (ReadPgfFile(fname & ".pgf") = True) Then

            If (ReadPulFile(fname & ".pul") = True) Then

                Dim b_reader As New BinaryReader(File.OpenRead(dir_name))
                Try
                    'READ .DAT FILE
                    ' for each group
                    org_xmin_value = 0
                    org_xmax_value = 0
                    For i = 0 To pul_rr.number_children - 1
                        ' for each series
                        For j = 0 To pul_rr.gr(i).number_children - 1
                            ' for each sweep
                            For k = 0 To pul_rr.gr(i).srr(j).number_children - 1
                                in_data(n + ns) = New InData(pul_rr.gr(i).srr(j).swr(k).total_points - 1)
                                in_data(n + ns).size = pul_rr.gr(i).srr(j).swr(k).total_points
                                in_data(n + ns).sweep = i + 1 & "_" & j + 1 & "_" & k + 1
                                in_data(n + ns).rate = pgf_rr.smr(j).sample_interval
                                in_data(n + ns).fname = dir_name
                                ' read the sound data first
                                For l = 0 To pul_rr.gr(i).srr(j).swr(k).total_points - 1
                                    tmp = b_reader.ReadInt16
                                    in_data(n + ns).Copy(Math.Round(l * in_data(n + ns).rate * 1000, 2), tmp / 1000000, tmp, l)
                                    'find(min And max)
                                    If (in_data(n + ns).id(l).trace1 > org_ymax_value) Then
                                        org_ymax_value = in_data(n + ns).id(l).trace1
                                    End If
                                    

                                    If (in_data(n + ns).id(l).trace1 < org_ymin_value) Then
                                        org_ymin_value = in_data(n + ns).id(l).trace1
                                    End If
                                Next
                                ' read the data
                                For l = 0 To pul_rr.gr(i).srr(j).swr(k).total_points - 1
                                    tmp = b_reader.ReadInt16
                                    tmp = tmp * pul_rr.gr(i).srr(j).swr(k).data_factor2
                                    in_data(n + ns).id(l).trace2 = tmp
                                    'find(min And max)
                                    If (in_data(n + ns).id(l).trace2 > org_ymax_value) Then
                                        org_ymax_value = in_data(n + ns).id(l).trace2
                                    End If

                                    If (in_data(n + ns).id(l).trace2 < org_ymin_value) Then
                                        org_ymin_value = in_data(n + ns).id(l).trace2
                                    End If
                                Next
                                If (in_data(n + ns).id(in_data(n + ns).size - 1).time > org_xmax_value) Then
                                    org_xmax_value = in_data(n + ns).id(in_data(n + ns).size - 1).time
                                End If
                                '' read the sound data first
                                'For l = 0 To pul_rr.gr(i).srr(j).swr(k).total_points - 1
                                '    tmp = b_reader.ReadInt16
                                'Next
                                '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                ns = ns + 1
                            Next
                        Next
                    Next
                Catch 'e As Exception
                    ns = 0
                    MsgBox(fname & ".dat" & ChrW(10) & ChrW(10) & "The file wasn't found or it's corrupted.")
                End Try

                b_reader.Close()
                i = 0 ' for reasons unknown to me without this dumb statement b_reader is not closed 
            Else
                MsgBox("The corresponding  .pul file wasn't " & ChrW(10) & "found or it's corrupted.")
                ns = 0
            End If
        Else
            MsgBox("The corresponding  .pgf file wasn't " & ChrW(10) & "found or it's corrupted.")
            ns = 0
        End If



        'in_data(0) = New InData(19)
        'in_data(0).size = 20
        'in_data(0).sweep = "1_1_1"
        'in_data(0).rate = 0.00001
        'in_data(0).fname = dir_name
        'in_data(0).Copy(0, -0.5, 0)
        'in_data(0).Copy(0.1, -0.4, 1)
        'in_data(0).Copy(0.28, -0.38, 2)
        'in_data(0).Copy(0.3, -0.3, 3)
        'in_data(0).Copy(0.4, -0.6, 4)
        'in_data(0).Copy(0.5, -0.3, 5)
        'in_data(0).Copy(0.6, -0.28, 6)
        'in_data(0).Copy(0.7, -0.4, 7)
        'in_data(0).Copy(0.8, -0.3, 8)
        'in_data(0).Copy(0.9, -0.26, 9)
        'in_data(0).Copy(1.0, -0.5, 10)
        'in_data(0).Copy(1.1, -0.55, 11)
        'in_data(0).Copy(1.2, -0.45, 12)
        'in_data(0).Copy(1.3, -0.35, 13)
        'in_data(0).Copy(1.4, -0.5, 14)
        'in_data(0).Copy(1.5, -0.3, 15)
        'in_data(0).Copy(1.6, -0.28, 16)
        'in_data(0).Copy(1.7, -0.5, 17)
        'in_data(0).Copy(1.8, -0.4, 18)
        'in_data(0).Copy(1.9, -0.3, 19)
        'org_ymin_value = -0.6
        'org_ymax_value = -0.26
        'org_xmin_value = 0
        'org_xmax_value = 1.9
        cur_xmin_value = org_xmin_value
        cur_xmax_value = org_xmax_value
        cur_ymin_value = org_ymin_value
        cur_ymax_value = org_ymax_value
        Return ns
    End Function

    Private Function ReadPgfFile(ByVal fn As String) As Boolean

        Dim i, j As Integer
        Dim rv As Boolean = True
        Dim size_of As Integer

        Using b_reader As New BinaryReader(File.OpenRead(fn))
            Try
                'READ .PGF FILE ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 
                ' scan tree
                pgf_tr.magic_number = b_reader.ReadInt32
                pgf_tr.number_tree_levels = b_reader.ReadInt32
                pgf_tr.RDim(pgf_tr.number_tree_levels - 1)
                For i = 0 To pgf_tr.number_tree_levels - 1
                    pgf_tr.tree_level(i) = b_reader.ReadInt32
                Next

                pgf_rr.root_version = b_reader.ReadInt16
                pgf_rr.number_children = b_reader.ReadInt32   ' number of  children - stimulations

                pgf_rr.RDim(pgf_rr.number_children - 1)
                ' read stimulation records
                For i = 0 To pgf_rr.number_children - 1
                    size_of = 0
                    pgf_rr.smr(i).file_name14 = ReadString(b_reader, STR_SIZE)
                    size_of += STR_SIZE
                    pgf_rr.smr(i).entry_name14 = ReadString(b_reader, STR_SIZE)
                    size_of += STR_SIZE
                    pgf_rr.smr(i).sample_interval = b_reader.ReadDouble ' sampling rate period
                    size_of += 8
                    pgf_rr.smr(i).filter_factor = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).st_sweep_interval = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).number_sweeps = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).number_repeats = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).repeat_wait = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).linked_sequence14 = ReadString(b_reader, STR_SIZE)
                    size_of += STR_SIZE
                    pgf_rr.smr(i).linked_wait = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).leak_count = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).leak_size = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).leak_holding = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).leak_alternate = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).alt_leak_averaging = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).leak_delay = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_segment1 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).trigger_time1 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_length1 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_amplitude1 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_dac1 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).trigger_segment2 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).trigger_time2 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_length2 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_amplitude2 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_dac2 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).trigger_segment3 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).trigger_time3 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_length3 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_amplitude3 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).trigger_dac3 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).number_of_triggers = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).relevant_x_segment = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).relevant_y_segment = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).write_mode = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).increment_mode = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).total_sweep_length = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).max_sweep_length = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).input_channels = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).g_update = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).rel_abs_pot = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).has_continuous = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).log_increment = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).pgf_stim_dac = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).adc1 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).adc2 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).y_unit1 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).y_unit2 = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).v_memb_increment = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).ext_trigger = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).file_template = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).stim_kind = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).lock_in_cycle = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).lock_in_amplitude = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_on = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).v_memb_mode = b_reader.ReadByte     'pgf_st_vmemb_mode
                    size_of += 1
                    pgf_rr.smr(i).fura_tot_length = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_delay = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_length1 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_length2 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_wave_length0 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_wave_length1 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_wave_length2 = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).fura_repeats = b_reader.ReadInt16
                    size_of += 2
                    pgf_rr.smr(i).lock_in_skip = b_reader.ReadInt32
                    size_of += 4
                    pgf_rr.smr(i).lock_in_v_reversal = b_reader.ReadDouble
                    size_of += 8
                    pgf_rr.smr(i).lock_in_mode = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).lock_in_show = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).config_macro16 = ReadString(b_reader, MACRO_NAME_SIZE)
                    size_of += MACRO_NAME_SIZE
                    pgf_rr.smr(i).end_macro16 = ReadString(b_reader, MACRO_NAME_SIZE)
                    size_of += MACRO_NAME_SIZE
                    pgf_rr.smr(i).ampl_mode_kind = b_reader.ReadByte
                    size_of += 1
                    pgf_rr.smr(i).filler1 = b_reader.ReadByte
                    size_of += 1
                    buffer = b_reader.ReadBytes(pgf_tr.tree_level(1) - size_of) '??????????????????????????????????????????????????

                    pgf_rr.smr(i).number_children = b_reader.ReadInt32   ' number of children - segments
                    ' read stimulation records
                    pgf_rr.smr(i).Rdim(pgf_rr.smr(i).number_children)
                    For j = 0 To pgf_rr.smr(i).number_children - 1
                        size_of = 0
                        pgf_rr.smr(i).sgr(j).se_class = b_reader.ReadByte
                        size_of += 1
                        pgf_rr.smr(i).sgr(j).is_holding = b_reader.ReadByte
                        size_of += 1
                        pgf_rr.smr(i).sgr(j).voltage = b_reader.ReadDouble
                        size_of += 8
                        pgf_rr.smr(i).sgr(j).duration = b_reader.ReadDouble
                        size_of += 8
                        pgf_rr.smr(i).sgr(j).delta_v_factor = b_reader.ReadDouble
                        size_of += 8
                        pgf_rr.smr(i).sgr(j).delta_v_increment = b_reader.ReadDouble
                        size_of += 8
                        pgf_rr.smr(i).sgr(j).delta_t_factor = b_reader.ReadDouble
                        size_of += 8
                        pgf_rr.smr(i).sgr(j).delta_t_increment = b_reader.ReadDouble
                        size_of += 8
                        buffer = b_reader.ReadBytes(pgf_tr.tree_level(2) - size_of)
                        pgf_rr.smr(i).sgr(j).number_children = b_reader.ReadInt32
                    Next
                Next
            Catch
                rv = False
            End Try
        End Using

        Return rv
    End Function

    Private Function ReadPulFile(ByVal fn As String) As Boolean

        Dim i, j, k As Integer
        Dim rv As Boolean = True
        Dim size_of As Integer

        Using b_reader As New BinaryReader(File.OpenRead(fn))
            Try
                'READ .PUL FILE '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                ' scan tree
                pul_tr.magic_number = b_reader.ReadInt32
                pul_tr.number_tree_levels = b_reader.ReadInt32

                pul_tr.RDim(pul_tr.number_tree_levels - 1)
                For i = 0 To pul_tr.number_tree_levels - 1
                    pul_tr.tree_level(i) = b_reader.ReadInt32
                Next

                'read root_record
                pul_rr.version = b_reader.ReadInt16
                pul_rr.version_name14 = ReadString(b_reader, STR_SIZE)
                pul_rr.file_name14 = ReadString(b_reader, STR_SIZE)
                pul_rr.comments400 = ReadString(b_reader, ROOT_TEXT_SIZE)
                pul_rr.start_time = b_reader.ReadDouble

                pul_rr.number_children = b_reader.ReadInt32   'number of children   -   groups
                pul_rr.Rdim(pul_rr.number_children - 1)

                ' read group records
                For i = 0 To pul_rr.number_children - 1
                    pul_rr.gr(i).label14 = ReadString(b_reader, STR_SIZE)
                    pul_rr.gr(i).text80 = ReadString(b_reader, COMMENT_SIZE)
                    pul_rr.gr(i).experimenter_number = b_reader.ReadInt32
                    pul_rr.gr(i).extra_long_real = b_reader.ReadDouble

                    pul_rr.gr(i).number_children = b_reader.ReadInt32     'number of children   -   series
                    pul_rr.gr(i).Rdim(pul_rr.gr(i).number_children - 1)

                    ' read series records
                    For j = 0 To pul_rr.gr(i).number_children - 1
                        size_of = 0
                        pul_rr.gr(i).srr(j).time = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).bandwidth = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).pipette_potential = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).cell_potential = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).pipette_resistance = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).seal_resistance = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).background_noise = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).temperature = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).pipette_pressure = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).user_param1_value = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).user_param1_name14 = ReadString(b_reader, STR_SIZE)
                        size_of += STR_SIZE
                        pul_rr.gr(i).srr(j).user_param1_unit = b_reader.ReadInt16
                        size_of += 2
                        pul_rr.gr(i).srr(j).user_param2_value = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).user_param2_name14 = ReadString(b_reader, STR_SIZE)
                        size_of += STR_SIZE
                        pul_rr.gr(i).srr(j).user_param2_unit = b_reader.ReadInt16
                        size_of += 2
                        pul_rr.gr(i).srr(j).recording_mode = b_reader.ReadByte
                        size_of += 1
                        pul_rr.gr(i).srr(j).filler1 = b_reader.ReadByte
                        size_of += 1
                        pul_rr.gr(i).srr(j).comment80 = ReadString(b_reader, COMMENT_SIZE)
                        size_of += COMMENT_SIZE
                        'pul_rr.gr(i).srr(j).epc9_state104 = ReadString(b_reader, EPC9_STATS_SIZE)
                        'size_of += EPC9_STATS_SIZE
                        pul_rr.gr(i).srr(j).internal_solution = b_reader.ReadInt32
                        size_of += 4
                        pul_rr.gr(i).srr(j).external_solution = b_reader.ReadInt32
                        size_of += 4
                        pul_rr.gr(i).srr(j).extra_y_unit1 = b_reader.ReadInt16
                        size_of += 2
                        pul_rr.gr(i).srr(j).extra_y_unit2 = b_reader.ReadInt16
                        size_of += 2
                        pul_rr.gr(i).srr(j).disp_y_unit1 = b_reader.ReadInt32
                        size_of += 4
                        pul_rr.gr(i).srr(j).disp_y_unit2 = b_reader.ReadInt32
                        size_of += 4
                        pul_rr.gr(i).srr(j).fura_k = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).fura_min = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).fura_max = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).lock_in_ext_phase = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).timer = b_reader.ReadDouble
                        size_of += 8
                        pul_rr.gr(i).srr(j).extra_long_real = b_reader.ReadDouble
                        size_of += 8
                        buffer = b_reader.ReadBytes(pul_tr.tree_level(2) - size_of)
                        pul_rr.gr(i).srr(j).number_children = b_reader.ReadInt32   'number of children   -   sweeps
                        pul_rr.gr(i).srr(j).Rdim(pul_rr.gr(i).srr(j).number_children)
                        ' read sweep records
                        For k = 0 To pul_rr.gr(i).srr(j).number_children - 1
                            size_of = 0
                            pul_rr.gr(i).srr(j).swr(k).time = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).stim_count = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).sweep_count = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).average_count = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).leak = b_reader.ReadByte
                            size_of += 1
                            pul_rr.gr(i).srr(j).swr(k).second_trace = b_reader.ReadByte
                            size_of += 1
                            pul_rr.gr(i).srr(j).swr(k).label14 = ReadString(b_reader, STR_SIZE)
                            size_of += STR_SIZE
                            pul_rr.gr(i).srr(j).swr(k).data_points = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).data = b_reader.ReadInt32
                            size_of += 4.0
                            pul_rr.gr(i).srr(j).swr(k).data_pointer = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).data_factor1 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).data_factor2 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).c_slow = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).g_series = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).rs_value = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).rs_fraction = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).zero_current = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).online_y_result = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).online_x_result = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).total_points = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).offset = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).sweep_kind = b_reader.ReadInt16
                            size_of += 2
                            pul_rr.gr(i).srr(j).swr(k).fura_points = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).fura_data = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).fura_pointer = b_reader.ReadInt32
                            size_of += 4
                            pul_rr.gr(i).srr(j).swr(k).online_y_result2 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).online_x_result2 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).disp_factor1 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).disp_factor2 = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).data_format = b_reader.ReadByte
                            size_of += 1
                            pul_rr.gr(i).srr(j).swr(k).data_abscissa = b_reader.ReadByte
                            size_of += 1
                            pul_rr.gr(i).srr(j).swr(k).timer = b_reader.ReadDouble
                            size_of += 8
                            pul_rr.gr(i).srr(j).swr(k).sweep_spares = b_reader.ReadInt32
                            size_of += 4
                            buffer = b_reader.ReadBytes(pul_tr.tree_level(3) - size_of)
                            pul_rr.gr(i).srr(j).swr(k).number_children = b_reader.ReadInt32       'number of children   -     = 0  now
                        Next
                    Next
                Next
            Catch
                rv = False
            End Try
        End Using

        Return rv
    End Function

    Private Shared Function ReadString(ByRef br As BinaryReader, ByVal size As Short) As String
        Dim char_array(size - 1) As Char
        Dim inta(size - 1) As Integer
        Dim str As String

        char_array = br.ReadChars(size)
        str = System.Convert.ToString(char_array)
        Return str
    End Function


End Class

