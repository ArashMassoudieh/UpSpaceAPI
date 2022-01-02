#include "Grid.h"
#include "StringOP.h"
#include "MapAsTimeSeriesSet.h"

void CGrid::runcommands_qt()
{

	set_progress_value(0);
	vector<vtkSmartPointer<vtkActor>> actors;
	CBTCSet All_Breakthroughpoints;
	CBTCSet Allpoints_velocities_eulerian(3);
	CBTCSet velocity_samples_from_trajs(3);
	CBTCSet extracted_OU_parameters;

	for (int i = 0; i < commands.size(); i++)
	{
            #ifdef QT_version
            main_window->setCursor(Qt::WaitCursor);
            QApplication::processEvents();
            qDebug() << QString::fromStdString(commands[i].command) << endl;
            #endif // QT_version
            cout<< "******* Running line #" << i << " *************" <<endl;

            if (commands[i].command == "get_all_btc_points")
            {
                show_in_window("Getting all btc points");
                All_Breakthroughpoints.getfromfile(pathout + commands[i].parameters["filename"],1);
            }

            if (commands[i].command == "get_points_v_eulerian")
            {
                show_in_window("Getting all Eulerian velocities");
                Allpoints_velocities_eulerian.getfromfile(pathout + commands[i].parameters["filename"],1);
            }
            if (commands[i].command == "get_points_v_from_traj")
            {
                show_in_window("Getting all Lagrangian velocities");
                velocity_samples_from_trajs.getfromfile(pathout + commands[i].parameters["filename"],1);
            }
            if (commands[i].command == "get_ou_params")
            {
                show_in_window("Getting all OU parameters");
                extracted_OU_parameters.getfromfile(pathout + commands[i].parameters["filename"],1);
            }


            if (commands[i].command == "write_all_btc_points")
            {
                show_in_window("Writing all btc points");
                All_Breakthroughpoints.writetofile(pathout + commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_points_v_eulerian")
            {
                show_in_window("Writing all Eulerian velocities");
                Allpoints_velocities_eulerian.writetofile(pathout + commands[i].parameters["filename"]);
            }
            if (commands[i].command == "write_points_v_from_traj")
            {
                show_in_window("Writing all Lagrangian velocities");
                velocity_samples_from_trajs.writetofile(pathout + commands[i].parameters["filename"]);
            }
            if (commands[i].command == "write_ou_params")
            {
                 show_in_window("Writing all OU parameters");
                extracted_OU_parameters.writetofile(pathout + commands[i].parameters["filename"]);
            }


            if (commands[i].command == "generate_k_field")
            {
                clear_contents();
                show_in_window("Assigning K...,");

                //cout << "Assigning K..." << endl;
                field_gen.max_correl_n = atoi(commands[i].parameters["n_neighbors"].c_str());
                field_gen.k_correlation_lenght_scale_x = atof(commands[i].parameters["corr_length_scale_x"].c_str());
                field_gen.k_correlation_lenght_scale_y = atof(commands[i].parameters["corr_length_scale_y"].c_str());
                if (commands[i].parameters.count("layered"))
                    {
                        if (commands[i].parameters["layered"] == "1")
                        {

                        }
                        else
                        assign_K_gauss();
                    }
                else
                    assign_K_gauss();
            }

            if (commands[i].command == "assign_linear_velocity_field")
            {
                show_in_window("assigning linear velocity field...");
                double v0 = atof(commands[i].parameters["v0"].c_str());
                double slope = atof(commands[i].parameters["v_slope"].c_str());
                Assign_Linear_Velocity_Field(v0, slope);
            }

            if (commands[i].command == "write_k_field")
            {
                show_in_window("Writing K field...");
                cout << "Writing K field..." << endl;
                writeasmatrixK(pathout+commands[i].parameters["filename"], 0);
            }

            if (commands[i].command == "create_1d_grid")
            {
                show_in_window("Creating one-D grid...");
                onedgrid.Generate_omega_field(atoi(commands[i].parameters["nx"].c_str()),atof(commands[i].parameters["dx"].c_str()),atof(commands[i].parameters["correlation_ls"].c_str()));
                show_in_window("One-D grid created!");
            }

            if (commands[i].command == "create_1d_grid_u")
            {
                show_in_window("Creating one-D grid u...");
                onedgrid.Generate_u_field(atoi(commands[i].parameters["nx"].c_str()));
                show_in_window("One-D grid created!");
            }

            if (commands[i].command == "write_1d_grid")
            {
                show_in_window("Writing one-D grid...");
                onedgrid.write_omega_field(pathout + commands[i].parameters["filename"]);

            }

            if (commands[i].command == "write_omega_dist")
            {
                show_in_window("Writing one-D grid...");
                onedgrid.getDist(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["filename"]);

            }

            if (commands[i].command == "write_u_dist_1d")
            {
                show_in_window("Writing u distribution on a one-D grid...");
                onedgrid.getDistU(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["filename"]);

            }


            if (commands[i].command == "assign_concentration_to_1d_grid")
            {
                show_in_window("Assigning concentration to 1d grid...");
                if (commands[i].parameters.count("i")>0)
                {
                    onedgrid.AssignConcentration(atoi(commands[i].parameters["i"].c_str()),atof(commands[i].parameters["value"].c_str()));
                }
                else
                    onedgrid.AssignConcentration(atof(commands[i].parameters["value"].c_str()));

            }

            if (commands[i].command == "assign_concentration_based_on_rank")
            {
                show_in_window("Assigning concentration based on rank");
                    onedgrid.AssignConcentration_based_on_u(atof(commands[i].parameters["u_min"].c_str()),atof(commands[i].parameters["u_max"].c_str()),atof(commands[i].parameters["value"].c_str()));

            }

            if (commands[i].command == "assign_concentration_to_u_grid")
            {
                show_in_window("Assigning concentration to u grid...");
                onedgrid.Initialize_U(atoi(commands[i].parameters["i"].c_str()),atof(commands[i].parameters["value"].c_str()));
            }

            if (commands[i].command == "write_c_dist_over_omega")
            {
                show_in_window("Writing C-dist...");
                onedgrid.GetConcentrationDistributionOverOmega(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["filename"]);
            }

            if (commands[i].command == "solve_1d")
            {
                show_in_window("Solving over x ...");
                onedgrid.Solve(atof(commands[i].parameters["dt"].c_str()),atof(commands[i].parameters["t_end"].c_str()),atof(commands[i].parameters["diffusion"].c_str()),atoi(commands[i].parameters["nbins"].c_str()));
                if (commands[i].parameters.count("filename_cx")>0)
                    onedgrid.ANSCX.writetofile(pathout + commands[i].parameters["filename_cx"],1,true);

                if (commands[i].parameters.count("filename_c_omega")>0)
                    onedgrid.ANSCW.Transpose(atof(commands[i].parameters["dt"].c_str()),"X").writetofile(pathout + commands[i].parameters["filename_c_omega"],1,true);

                if (commands[i].parameters.count("filename_c_u")>0)
                    onedgrid.ANSCU.Transpose(atof(commands[i].parameters["dt"].c_str()),"X").writetofile(pathout + commands[i].parameters["filename_c_u"],1,true);

                if (commands[i].parameters.count("filename_u_dev")>0)
                    onedgrid.u_dev.writefile(pathout + commands[i].parameters["filename_u_dev"]);


                show_in_window("Solving over x, done!");
            }

            if (commands[i].command == "solve_1d_u")
            {
                show_in_window("Solving over u ...");
                onedgrid.Solve_U(atof(commands[i].parameters["dt"].c_str()),atof(commands[i].parameters["t_end"].c_str()),atof(commands[i].parameters["diffusion"].c_str()),atof(commands[i].parameters["correlation"].c_str()));
                if (commands[i].parameters.count("filename_c_u")>0)
                    onedgrid.ANSCU.writetofile(pathout + commands[i].parameters["filename_c_u"],1,true);

                 if (commands[i].parameters.count("filename_u_dev")>0)
                    onedgrid.u_dev.writefile(pathout + commands[i].parameters["filename_u_dev"]);
                //if (commands[i].parameters.count("filename_c_omega")>0)
                //    onedgrid.ANSCW.Transpose(atof(commands[i].parameters["dt"].c_str()),"X").writetofile(pathout + commands[i].parameters["filename_c_omega"],1,true);

                //if (commands[i].parameters.count("filename_c_u")>0)
                //    onedgrid.ANSCU.Transpose(atof(commands[i].parameters["dt"].c_str()),"X").writetofile(pathout + commands[i].parameters["filename_c_u"],1,true);


                show_in_window("Solving over u, done!");
            }



            if (commands[i].command == "read_k_field")
            {
                    show_in_window("Reading K field...");
                    cout << "Reading K field..." << endl;
                    GP.nx = atoi(commands[i].parameters["nx"].c_str());
                    GP.ny = atoi(commands[i].parameters["ny"].c_str());
                    GP.nx_data = atoi(commands[i].parameters["nx_data"].c_str());
                    GP.ny_data = atoi(commands[i].parameters["ny_data"].c_str());
                    getKfromfile(pathout+commands[i].parameters["filename"], GP.nx, GP.ny, GP.nx_data, GP.ny_data);
            }
            if (commands[i].command == "solve_hydro")
            {
                    show_in_window("Solving hydro ...");
                    leftboundary_h = atof(commands[i].parameters["l_boundary"].c_str());
                    rightboundary_h = atof(commands[i].parameters["r_boundary"].c_str());
                    //cout << "Solving hydro ..." << endl;
                    CMatrix H = solve();
                    Assign_Ranks();
            }

            if (commands[i].command == "solve_transport")
            {
                    show_in_window("Solving transport ...");
                    numberofspecies = 1;
                    if (commands[i].parameters.count("nspecies")>0)
                    {
                        numberofspecies = atoi(commands[i].parameters["nspecies"].c_str());
                    }
                    time_weight = atof(commands[i].parameters["weight"].c_str());
                    leftboundary_C = ATOF(split(commands[i].parameters["l_boundary"],','));
                    if (leftboundary_C.size()==0) leftboundary_C.push_back(1);
                    D = atof(commands[i].parameters["diffusion"].c_str());
                    dt = atof(commands[i].parameters["dt"].c_str());

                    vector<double> decay_coeff(numberofspecies);
                    vector<double> decay_order(numberofspecies);
                    if (commands[i].parameters.count("decay_coeff")>0)
                        decay_coeff = ATOF(split(commands[i].parameters["decay_coeff"],','));

                    if (commands[i].parameters.count("decay_order")>0)
                        decay_order = ATOF(split(commands[i].parameters["decay_order"],','));

                    cout << "Solving transport ..." << endl;
                    solve_transport(atof(commands[i].parameters["t_end"].c_str()),decay_coeff,decay_order);
            }

            if (commands[i].command == "write_btc_from_concentration")
            {
                int species_id=0;
                if (commands[i].parameters.count("species")>0)
                    species_id = atoi(commands[i].parameters["species"].c_str());
                show_in_window("Writing Breakthrough curve at x = " + commands[i].parameters["x"] + "...");
                GetConcentrationBTCAtX(species_id, atof(commands[i].parameters["x"].c_str()),commands[i].parameters["filename"],commands[i].parameters["filename_d"]);

            }

            if (commands[i].command == "solve_transport_ou")
            {
                show_in_window("Solving transport (Ornstein-Uhlenbeck)...");
                time_weight = atof(commands[i].parameters["weight"].c_str());
                D = atof(commands[i].parameters["diffusion"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());
                OU.lc = atof(commands[i].parameters["lc"].c_str());
                OU.ld = atof(commands[i].parameters["ld"].c_str());
                OU.diffusion = atof(commands[i].parameters["diffusion"].c_str());
                double decay_coeff = atof(commands[i].parameters["decay_coeff"].c_str());
                double decay_order = atof(commands[i].parameters["decay_order"].c_str());
                solve_transport_OU(atof(commands[i].parameters["t_end"].c_str()),decay_coeff, decay_order);
                if (commands[i].parameters.count("filename") > 0) OU.BTCs.writetofile(pathout + commands[i].parameters["filename"]);
                if (commands[i].parameters.count("filename_d") > 0) OU.BTCs.detivative().writetofile(pathout + commands[i].parameters["filename_d"]);
                if (commands[i].parameters.count("filename_n") > 0) OU.BTC_normal.writetofile(pathout + commands[i].parameters["filename_n"]);
                if (commands[i].parameters.count("filename_nd") > 0) OU.BTC_normal.detivative().writetofile(pathout + commands[i].parameters["filename_nd"]);
                if (commands[i].parameters.count("filename_fw_nd") > 0) OU.BTC_normal_fw.detivative().writetofile(pathout + commands[i].parameters["filename_fw_nd"]);
            }

            if (commands[i].command == "solve_transport_copula")
            {
                show_in_window("Solving transport (Copula)...");
                time_weight = atof(commands[i].parameters["weight"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());
                numberofspecies = 1;
                if (commands[i].parameters.count("nspecies")>0)
                {
                    numberofspecies = atoi(commands[i].parameters["nspecies"].c_str());
                }
                time_weight = atof(commands[i].parameters["weight"].c_str());
                leftboundary_C = ATOF(split(commands[i].parameters["l_boundary"],','));
                if (leftboundary_C.size()==0) leftboundary_C.push_back(1);
                double diffusion = atof(commands[i].parameters["diffusion"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());

                vector<double> decay_coeff(numberofspecies);
                vector<double> decay_order(numberofspecies);
                if (commands[i].parameters.count("decay_coeff")>0)
                    decay_coeff = ATOF(split(commands[i].parameters["decay_coeff"],','));

                if (commands[i].parameters.count("decay_order")>0)
                    decay_order = ATOF(split(commands[i].parameters["decay_order"],','));
                copula_params.epsilon = atof(commands[i].parameters["epsilon"].c_str());
                copula_params.tau = atof(commands[i].parameters["tau"].c_str());
                copula_params.diffusion = atof(commands[i].parameters["diffusion"].c_str());
                copula_params.mean_method = commands[i].parameters["mean_method"];
                solve_transport_Copula(atof(commands[i].parameters["t_end"].c_str()),diffusion,decay_coeff, decay_order);
                if (commands[i].parameters.count("filename") > 0) OU.BTCs.writetofile(pathout + commands[i].parameters["filename"]);
                if (commands[i].parameters.count("filename_d") > 0) OU.BTCs.detivative().writetofile(pathout + commands[i].parameters["filename_d"]);
                if (commands[i].parameters.count("filename_n") > 0) OU.BTC_normal.writetofile(pathout + commands[i].parameters["filename_n"]);
                if (commands[i].parameters.count("filename_nd") > 0) OU.BTC_normal.detivative().writetofile(pathout + commands[i].parameters["filename_nd"]);
                if (commands[i].parameters.count("filename_fw_nd") > 0) OU.BTC_normal_fw.detivative().writetofile(pathout + commands[i].parameters["filename_fw_nd"]);
            }

            if (commands[i].command == "solve_transport_copula_diffusion")
            {
                show_in_window("Solving transport (Copula)...");
                time_weight = atof(commands[i].parameters["weight"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());
                numberofspecies = 1;
                if (commands[i].parameters.count("nspecies")>0)
                {
                    numberofspecies = atoi(commands[i].parameters["nspecies"].c_str());
                }
                time_weight = atof(commands[i].parameters["weight"].c_str());
                leftboundary_C = ATOF(split(commands[i].parameters["l_boundary"],','));
                double diffusion = atof(commands[i].parameters["diffusion"].c_str());
                dt = atof(commands[i].parameters["dt"].c_str());

                vector<double> decay_coeff(numberofspecies);
                vector<double> decay_order(numberofspecies);
                if (commands[i].parameters.count("decay_coeff")>0)
                    decay_coeff = ATOF(split(commands[i].parameters["decay_coeff"],','));

                if (commands[i].parameters.count("decay_order")>0)
                    decay_order = ATOF(split(commands[i].parameters["decay_order"],','));
                copula_params.epsilon = atof(commands[i].parameters["epsilon"].c_str());
                copula_params.tau = atof(commands[i].parameters["tau"].c_str());
                copula_params.diffusion = atof(commands[i].parameters["diffusion"].c_str());
                copula_params.mean_method = commands[i].parameters["mean_method"];
                solve_transport_Copula_diffusion(atof(commands[i].parameters["t_end"].c_str()),diffusion,decay_coeff, decay_order);
                if (commands[i].parameters.count("filename") > 0) OU.BTCs.writetofile(pathout + commands[i].parameters["filename"]);
                if (commands[i].parameters.count("filename_d") > 0) OU.BTCs.detivative().writetofile(pathout + commands[i].parameters["filename_d"]);
                if (commands[i].parameters.count("filename_n") > 0) OU.BTC_normal.writetofile(pathout + commands[i].parameters["filename_n"]);
                if (commands[i].parameters.count("filename_nd") > 0) OU.BTC_normal.detivative().writetofile(pathout + commands[i].parameters["filename_nd"]);
                if (commands[i].parameters.count("filename_fw_nd") > 0) OU.BTC_normal_fw.detivative().writetofile(pathout + commands[i].parameters["filename_fw_nd"]);
            }


            if (commands[i].command == "solve_transport_laplace")
            {
                show_in_window("Solving transport (laplace)...");
                for (int ii=0; ii<numberofspecies; ii++)
                    leftboundary_C[ii] = 1;
                D = atof(commands[i].parameters["diffusion"].c_str());

                solve_transport_laplace(atof(commands[i].parameters["s"].c_str()));
            }

            if (commands[i].command == "write_h_field")
            {
                show_in_window("Writing H field");
               #ifdef QT_version
                QApplication::processEvents();
               #endif // QT_version
                H.writetofile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_v_field")
            {
                show_in_window("Writing velocities ...");

                //cout << "Writing velocities ..." << endl;
                writeasmatrix(pathout+commands[i].parameters["filename_x"], 0);
                writeasmatrix(pathout+commands[i].parameters["filename_y"], 1);
            }

            if (commands[i].command == "initialize_trajectories")
            {
                show_in_window("Initializing trajectories ...");
                cout << "Initializing trajectories ..." << endl;
                bool flow_weighted;
                if (commands[i].parameters.count("flow_weighted")>0)
                    flow_weighted = atoi(commands[i].parameters["flow_weighted"].c_str());
                else
                    flow_weighted = true;

                CBTC ini_vel;
                if (commands[i].parameters["boundary_v_dist_filename"]!="")
                    ini_vel=initialize(atoi(commands[i].parameters["n"].c_str()), atof(commands[i].parameters["x_0"].c_str()), 0, flow_weighted, pathout+commands[i].parameters["boundary_v_dist_filename"], atoi(commands[i].parameters["weighted"].c_str()));
                else
                    ini_vel=initialize(atoi(commands[i].parameters["n"].c_str()), atof(commands[i].parameters["x_0"].c_str()), 0, flow_weighted, "", atoi(commands[i].parameters["weighted"].c_str()));
                weighted = atoi(commands[i].parameters["weighted"].c_str());
                cout << "weighted = " << weighted;
                if (commands[i].parameters.count("filename"))
                    ini_vel.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "create_trajectories")
            {
                show_in_window("Simulating trajectories ...");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                cout << "Simulating trajectories ..." << endl;
                Traj = gettrajectories_vdt(atof(commands[i].parameters["dt"].c_str()), atof(commands[i].parameters["x_end"].c_str()),atof(commands[i].parameters["tol"].c_str()),atof(commands[i].parameters["diffusion"].c_str()));
            }

            if (commands[i].command == "create_trajectories_fix_dx")
            {
                show_in_window("Simulating trajectories with fixed dx...");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                cout << "Simulating trajectories with fixed dx..." << endl;
                Traj = gettrajectories_fixed_dx(atof(commands[i].parameters["dx"].c_str()), atof(commands[i].parameters["x_end"].c_str()), atof(commands[i].parameters["tol"].c_str()), atof(commands[i].parameters["diffusion"].c_str()));
            }

            if (commands[i].command == "write_breakthrough_curve")
            {
                show_in_window("Get breakthrough curve at x = " + commands[i].parameters["x"]);
                CBTC Breakthroughcurve_from_trajs = Traj.get_BTC(atof(commands[i].parameters["x"].c_str()), atoi(commands[i].parameters["nbins"].c_str()), atoi(commands[i].parameters["velweight"].c_str()), atof(commands[i].parameters["smoothing_factor"].c_str()));
                if (All_Breakthroughpoints.lookup("x=" + commands[i].parameters["x"]) == -1)
                {
                    CBTC BTC = Traj.get_BTC_points(atof(commands[i].parameters["x"].c_str()),atoi(commands[i].parameters["velweight"].c_str()));
                    All_Breakthroughpoints.append(BTC, "x=" + commands[i].parameters["x"]);

                }
                else
                {   CBTC BTC = Traj.get_BTC_points(atof(commands[i].parameters["x"].c_str()),atoi(commands[i].parameters["velweight"].c_str()));
                    All_Breakthroughpoints.BTC[All_Breakthroughpoints.lookup("x=" + commands[i].parameters["x"])].append(BTC);
                }
                Breakthroughcurve_from_trajs.getcummulative().writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "get_profile")
            {
                show_in_window("Getting profile ... ");
                string filename = pathout + commands[i].parameters["filename"];
                double x_start = atof(commands[i].parameters["x_start"].c_str());
                double x_end = atof(commands[i].parameters["x_end"].c_str());
                double interval = atof(commands[i].parameters["interval"].c_str());
                int species_id = 0;
                if (commands[i].parameters.count("species")>0)
                    species_id = atoi(commands[i].parameters["species"].c_str());
                int timestep;
                if (commands[i].parameters.count("at_time")==0)
                    timestep = p[0][0].C.size()-1;
                else
                    timestep = min(int(p[0][0].C.size()-1),int(atof(commands[i].parameters["at_time"].c_str())/dt));

                GetProfile(species_id, timestep, x_start, x_end, interval, filename);
            }

            if (commands[i].command == "write_breakthrough_curves_all")
            {
                show_in_window("Writing break through curves for all realizations");
                All_Breakthroughpoints.distribution(atoi(commands[i].parameters["nbins"].c_str()), All_Breakthroughpoints.BTC.size(),0, atof(commands[i].parameters["smoothing_factor"].c_str())).getcummulative().writetofile(pathout + commands[i].parameters["filename"],true);
            }



            if (commands[i].command == "write_trajectories")
            {
                show_in_window("Writing trajectories ...");
                //cout << "Writing trajectories ..." << endl;
                if (commands[i].parameters.count("interval")>0)
                    Traj.write(pathout+commands[i].parameters["filename"]);
                else
                    Traj.write(pathout+commands[i].parameters["filename"],atoi(commands[i].parameters["interval"].c_str()));
            }

            if (commands[i].command == "get_velocity_dist")
            {
                show_in_window("Get Velocities into vectors...");

                CBTC vx = get_v_btc(0);
                CBTC vy = get_v_btc(1);
                CBTC v_mag = get_v_mag_btc();
                Allpoints_velocities_eulerian.BTC[0].append(vx);
                Allpoints_velocities_eulerian.BTC[1].append(vy);
                Allpoints_velocities_eulerian.BTC[2].append(v_mag);

                show_in_window("Get Velocities distributions");
                #ifdef QT_version
                QApplication::processEvents();
                #endif // QT_version
                //cout << "Get Velocities distributions" << endl;
                vx_dist = vx.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vx.maxC()-vx.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vy_dist = vy.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vy.maxC()-vy.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                v_dist = v_mag.distribution(atoi(commands[i].parameters["nbins"].c_str()), (v_mag.maxC()-v_mag.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
            }

            if (commands[i].command == "write_velocity_dist")
            {
                show_in_window("Writing velocities distributions...");

                //cout << "Writing velocities into vectors" << endl;
                vx_dist.writefile(pathout+commands[i].parameters["filename_x"]);
                vy_dist.writefile(pathout+commands[i].parameters["filename_y"]);
                v_dist.writefile(pathout+commands[i].parameters["filename_mag"]);
            }

            if (commands[i].command == "write_velocity_dist_mf")
            {
                show_in_window("Reading spatial velocity...");
                CBTC vx;
                if (commands[i].parameters.count("v_filename")==1)
                    vx = get_v_dist_MODFlow(commands[i].parameters["v_filename"]);
                else
                    vx = Traj.sample_velocities();
                vx.distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["filename"]);
                CVector stats(2);
                stats[0] = vx.Log(1e-6).mean();
                stats[1] = vx.Log(1e-6).std();
                if (commands[i].parameters.count("stat_filename"))
                    stats.writetofile(pathout+commands[i].parameters["stat_filename"]);
            }

            if (commands[i].command == "write_velocity_dist_frac")
            {
                show_in_window("Reading spatial velocity...");

                CBTC  vx;
                if (commands[i].parameters.count("v_filename")==1)
                    vx = get_v_dist_frac(commands[i].parameters["v_filename"]);
                else
                    vx = Traj.sample_velocities();

                vx.distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["filename"]);
                CVector stats(2);
                stats[0] = vx.Log(1e-6).mean();
                stats[1] = vx.Log(1e-6).std();
                if (commands[i].parameters.count("stat_filename"))
                    stats.writetofile(pathout+commands[i].parameters["stat_filename"]);
            }

            if (commands[i].command == "get_btc_frac")
            {

                show_in_window("getting btc at " + commands[i].parameters["x_min"] + "-" + commands[i].parameters["x_max"]);
                CBTCSet btc;
                CBTCSet btc_react;
                bool consider_reactions=false;
                if (commands[i].parameters.count("consider_reactions")>0)
                {
                    if (tolower(commands[i].parameters["consider_reactions"])=="true")
                        consider_reactions = true;
                    else
                        consider_reactions = false;
                }
                if (commands[i].parameters.count("filename")==1)
                    btc = get_BTC_frac(commands[i].parameters["filename"],atof(commands[i].parameters["x_min"].c_str()),atof(commands[i].parameters["x_max"].c_str()));
                else
                {
                    btc = get_BTC_frac(Traj,atof(commands[i].parameters["x_min"].c_str()),atof(commands[i].parameters["x_max"].c_str()));
                    if (consider_reactions)
                        btc_react = get_BTC_frac(Traj,atof(commands[i].parameters["x_min"].c_str()),atof(commands[i].parameters["x_max"].c_str()),consider_reactions);
                    else
                        btc_react = btc;

                }
                double coefficient=1;
                if (consider_reactions)
                {
                    coefficient = double(btc_react.BTC[0].n)/double(btc.BTC[0].n);
                }
                if (commands[i].parameters.count("raw_data_filename")>0)
                {
                    show_in_window("writing raw data...");
                    btc.writetofile(pathout + commands[i].parameters["raw_data_filename"]);
                }

                if (commands[i].parameters.count("btc_filename")>0)
                {
                    show_in_window("writing btc ...");
                    if (commands[i].parameters["log"]=="1")
                        (btc_react.BTC[0].distribution_log(atoi(commands[i].parameters["nbins"].c_str()))*coefficient).writefile(pathout + commands[i].parameters["btc_filename"]);
                    else
                        (btc_react.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()))*coefficient).writefile(pathout + commands[i].parameters["btc_filename"]);

                    if (commands[i].parameters.count("cdf_filename")>0)
                    {
                        show_in_window("writing btc ...");
                        if (commands[i].parameters["log"]=="1")
                            (btc_react.BTC[0].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),true)*coefficient).writefile(pathout + commands[i].parameters["cdf_filename"]);
                        else
                            (btc_react.BTC[0].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),false)*coefficient).writefile(pathout + commands[i].parameters["cdf_filename"]);
                    }
                }

                if (commands[i].parameters.count("filename")==0)
                if (commands[i].parameters.count("v_filename")>0)
                {
                    show_in_window("writing v_dist ...");
                    btc.BTC[1].distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["v_filename"]);
                }

                if (commands[i].parameters.count("stat_filename"))
                {
                    show_in_window("writing correlation");
                    CVector stats(1);
                    stats[0] = R(btc.BTC[1],1.0/btc.BTC[0],0);
                    stats.writetofile(pathout+commands[i].parameters["stat_filename"]);
                }

            }

             if (commands[i].command == "get_btc_mf")
            {

                show_in_window("getting btc (mf) at " + commands[i].parameters["x_min"] + "-" + commands[i].parameters["x_max"]);
                CBTCSet btc;
                if (commands[i].parameters.count("filename")==1)
                    btc = get_BTC_mf(commands[i].parameters["filename"],atof(commands[i].parameters["x_min"].c_str()),atof(commands[i].parameters["x_max"].c_str()));
                else
                    btc = get_BTC_frac(Traj,atof(commands[i].parameters["x_min"].c_str()),atof(commands[i].parameters["x_max"].c_str()));

                if (commands[i].parameters.count("raw_data_filename")>0)
                {
                    show_in_window("writing raw data...");
                    btc.writetofile(pathout + commands[i].parameters["raw_data_filename"]);
                }

                if (commands[i].parameters.count("btc_filename")>0)
                {
                    show_in_window("writing btc ...");
                    if (commands[i].parameters["log"]=="1")
                        btc.BTC[0].distribution_log(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["btc_filename"]);
                    else
                        btc.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["btc_filename"]);

                    if (commands[i].parameters.count("cdf_filename")>0)
                    {
                        show_in_window("writing btc ...");
                        if (commands[i].parameters["log"]=="1")
                            btc.BTC[0].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),true).writefile(pathout + commands[i].parameters["cdf_filename"]);
                        else
                            btc.BTC[0].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),false).writefile(pathout + commands[i].parameters["cdf_filename"]);
                    }
                }

                if (commands[i].parameters.count("filename")==0)
                if (commands[i].parameters.count("v_filename")>0)
                {
                    show_in_window("writing v_dist ...");
                    btc.BTC[1].distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout + commands[i].parameters["v_filename"]);
                }

                if (commands[i].parameters.count("stat_filename"))
                {
                    show_in_window("writing correlation");
                    CVector stats(1);
                    stats[0] = R(btc.BTC[1],1.0/btc.BTC[0],0);
                    stats.writetofile(pathout+commands[i].parameters["stat_filename"]);
                }

            }


            if (commands[i].command == "write_velocity_dist_all")
            {
                show_in_window("Writing velocities distributions for all realizations...");

                CBTC vx_dist_all = Allpoints_velocities_eulerian.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].maxC()-Allpoints_velocities_eulerian.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC vy_dist_all = Allpoints_velocities_eulerian.BTC[1].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[1].maxC()-Allpoints_velocities_eulerian.BTC[1].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC v_mag_dist_all = Allpoints_velocities_eulerian.BTC[2].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].maxC()-Allpoints_velocities_eulerian.BTC[2].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));

                CVector stats(2);
                stats[0] = Allpoints_velocities_eulerian.BTC[0].Log(1e-6).mean();
                stats[1] = Allpoints_velocities_eulerian.BTC[0].Log(1e-6).std();
                if (commands[i].parameters.count("stat_filename"))
                    stats.writetofile(pathout+commands[i].parameters["stat_filename"]);

                vx_dist_all.writefile(pathout + commands[i].parameters["filename_x"]);
                vy_dist_all.writefile(pathout + commands[i].parameters["filename_y"]);
                v_mag_dist_all.writefile(pathout + commands[i].parameters["filename_mag"]);
                if (commands[i].parameters.count("filename_x_log") > 0)
                {
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions" << endl;
                    #else
                    cout << "Calculating log-transformed distributions" << endl;
                    #endif // QT_version
                    CBTC vx_dist_all_log = Allpoints_velocities_eulerian.BTC[0].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].Log().maxC()-Allpoints_velocities_eulerian.BTC[0].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions, Done!" << endl;
                    #else
                    cout<< "Calculating log-transformed distributions, Done!" << endl;
                    #endif // QT_version
                    vx_dist_all_log.writefile(pathout + commands[i].parameters["filename_x_log"]);

                    if (commands[i].parameters.count("filename_log_int") > 0)
                    {
                        vx_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int"]);
                    }
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions for velocity magnitude" << endl;
                    #else
                    cout<< "Calculating log-transformed distributions for velocity magnitude" << endl;
                    #endif // QT_version
                    CBTC v_mag_dist_all_log = Allpoints_velocities_eulerian.BTC[2].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].Log().maxC()-Allpoints_velocities_eulerian.BTC[2].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                    #ifdef QT_version
                    qDebug() << "Calculating log-transformed distributions for velocity magnitude, Done!" << endl;
                    #else
                    cout << "Calculating log-transformed distributions for velocity magnitude, Done!" << endl;
                    #endif // QT_version
                    v_mag_dist_all_log.writefile(pathout + commands[i].parameters["filename_mag_log"]);
                    if (commands[i].parameters.count("filename_log_int_mag") > 0)
                    {
                        v_mag_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int_mag"]);
                    }
                }
            }

            if (commands[i].command == "write_velocity_dist_all_fw")
            {
                show_in_window("Writing velocities distributions for all realizations (flux weighted)...");

                CBTC vx_dist_all = Allpoints_velocities_eulerian.BTC[0].distribution_fw(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].maxC()-Allpoints_velocities_eulerian.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC vy_dist_all = Allpoints_velocities_eulerian.BTC[1].distribution_fw(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[1].maxC()-Allpoints_velocities_eulerian.BTC[1].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                CBTC v_mag_dist_all = Allpoints_velocities_eulerian.BTC[2].distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].maxC()-Allpoints_velocities_eulerian.BTC[2].minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vx_dist_all.writefile(pathout + commands[i].parameters["filename_x"]);
                vy_dist_all.writefile(pathout + commands[i].parameters["filename_y"]);
                v_mag_dist_all.writefile(pathout + commands[i].parameters["filename_mag"]);
                if (commands[i].parameters.count("filename_x_log") > 0)
                {

                    CBTC vx_dist_all_log = Allpoints_velocities_eulerian.BTC[0].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[0].Log().maxC()-Allpoints_velocities_eulerian.BTC[0].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).make_flux_weighted("log");
                    vx_dist_all_log.writefile(pathout + commands[i].parameters["filename_x_log"]);
                    if (commands[i].parameters.count("filename_log_int") > 0)
                    {
                        vx_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int"]);
                    }

                    CBTC v_mag_dist_all_log = Allpoints_velocities_eulerian.BTC[2].Log().distribution(atoi(commands[i].parameters["nbins"].c_str()), (Allpoints_velocities_eulerian.BTC[2].Log().maxC()-Allpoints_velocities_eulerian.BTC[2].Log().minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).make_flux_weighted("log");
                    v_mag_dist_all_log.writefile(pathout + commands[i].parameters["filename_mag_log"]);
                    if (commands[i].parameters.count("filename_log_int_mag") > 0)
                    {
                        v_mag_dist_all_log.unlog().writefile(pathout + commands[i].parameters["filename_log_int_mag"]);
                    }
                }
            }

            if (commands[i].command == "get_velocity_dist_at_sects")
            {
                show_in_window("Get Velocities distributions at cross-sections ...");

                //cout << "Get Velocities distributions at cross-sections ..." << endl;

                CBTC vx = get_v_btc(atof(commands[i].parameters["x"].c_str()), 0);
                CBTC vy = get_v_btc(atof(commands[i].parameters["x"].c_str()), 1);
                vx_dist = vx.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vx.maxC()-vx.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                vy_dist = vy.distribution(atoi(commands[i].parameters["nbins"].c_str()), (vy.maxC()-vy.minC())*atof(commands[i].parameters["smoothing_factor"].c_str()));
                sect_dist.append(vx_dist, "x=" + commands[i].parameters["x"]);
                sect_dist.append(vy_dist, "x=" + commands[i].parameters["x"]);

            }

            if (commands[i].command == "write_velocities_at_sects")
            {
                show_in_window("Write velocities distributions at cross-sections ...");

                //cout << "Write velocities distributions at cross-sections ..." << endl;
                sect_dist.writetofile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "create_inverse_marginal_k")
            {
                show_in_window("Create inverse marginal hydraulic conductivity ...");

                //cout << "Create inverse marginal hydraulic conductivity ..." << endl;
                set_inv_K_dist(atoi(commands[i].parameters["ninc"].c_str()));
            }

            if (commands[i].command == "write_inverse_marginal_k")
            {
                show_in_window("Write inverse marginal hydraulic conductivity ...");

                //cout << "Write inverse marginal hydraulic conductivity ..." << endl;
                inv_K_dist.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_marginal_k")
            {
                show_in_window("Write marginal hydraulic conductivity ...");

                //cout << "Write marginal hydraulic conductivity ..." << endl;
                CBTC K = get_K_CDF(atof(commands[i].parameters["x0"].c_str()), atof(commands[i].parameters["x1"].c_str()), atof(commands[i].parameters["log_inc"].c_str()));
                K.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_marginal_v_dist")
            {
                show_in_window("Write marginal velocity dist ...");

                //cout << "Write marginal hydraulic conductivity ..." << endl;
                double x0 = atof(commands[i].parameters["x0"].c_str());
                bool _log = true;
                if (x0 < 0) _log = false;
                CBTC K = get_V_PDF(atof(commands[i].parameters["x0"].c_str()), atof(commands[i].parameters["x1"].c_str()), atof(commands[i].parameters["log_inc"].c_str()),_log);
                K.writefile(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "renormalize_k")
            {
                show_in_window("Renormalizing hydraulic conductivity ...");

                //cout << "Renormalizing hydraulic conductivity ..." << endl;
                renormalize_k();
            }

            #ifdef QT_version
            if (commands[i].command == "show_k_field")
            {
                show_in_window("Showing K field...");
                show_K_field();
            }
            #endif // QT_version

            if (commands[i].command == "save_k_field_into_vtp")
            {
                show_in_window("Saving K field (VTK)...");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_K_field_to_vtp(pathout+commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()),atoi(commands[i].parameters["offset"].c_str()));
            }

            if (commands[i].command == "save_solution_into_vtp")
            {

                show_in_window("Saving solution (VTK)...");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_K_solution_to_vtp(pathout + commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atoi(commands[i].parameters["log"].c_str()));
            }

            if (commands[i].command == "save_concentration_into_vtp")
            {
                int species_counter;
                if (numberofspecies==1)
                    species_counter = 0;
                else if (commands[i].parameters.count("species")==0)
                    species_counter = 0;
                else
                    species_counter = atoi(commands[i].parameters["species"].c_str());
                double concentration_interval = -1;
                show_in_window("Saving Concentration (VTK)...");
                if (commands[i].parameters.count("interval") > 0)
                    concentration_interval = atof(commands[i].parameters["interval"].c_str());
                else
                    concentration_interval = 1;
                vector<double> intervals;
                for (int i = 0; i < p[0][0].C.size(); i+=concentration_interval)
                {
                        intervals.push_back(i*dt);
                }
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "surface.vtp";
                write_C_to_vtp(pathout + commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atoi(commands[i].parameters["log"].c_str()), intervals);

            }

            if (commands[i].command == "show")
            {
                show_in_window("Showing so far...");

                if (commands[i].parameters.count("filename") == 0)
                        showthings(actors);
                else if (commands[i].parameters.size() == 1)
                        showthings(actors,pathout+commands[i].parameters["filename"]);
                actors.clear();
            }

            if (commands[i].command == "get_trajectory")
            {
                show_in_window("Showing trajectory " + commands[i].parameters[0]);
                actors.push_back(traj_vtk_pdt(atoi(commands[i].parameters["traj_no"].c_str()), atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str())));
            }

            if (commands[i].command == "get_trajectories")
            {
                show_in_window("Converting trajectories into graphic objects... ");
                vector<vtkSmartPointer<vtkActor>> tmp_actor;
                tmp_actor = trajs_vtk_pdt(atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()));

                for (int i = 0; i < tmp_actor.size(); i++)
                        actors.push_back(tmp_actor[i]);

            }

            #if QT_version
            if (commands[i].command == "screen_shot_test")
            {
                show_in_window("Screen Shot test");
                screenshot_test();
            }
            #endif // QT_version

            if (commands[i].command == "save_trajs_into_vtp")
            {
                show_in_window("Writing trajectories into vtp... ");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "paths.vtp";
                int interval = 1;
                if (commands[i].parameters.count("interval") == 0)
                    interval = 1;
                else
                    interval = atoi(commands[i].parameters["interval"].c_str());
                trajs_vtk_pdt_to_vtp(pathout+commands[i].parameters["filename"], atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()), atof(commands[i].parameters["log"].c_str()), atof(commands[i].parameters["color"].c_str()),interval);
            }

            if (commands[i].command == "save_trajs_into_vtp_3d")
            {
                show_in_window("Writing trajectories into vtp... ");
                if (commands[i].parameters.count("filename") == 0) commands[i].parameters["filename"] == "paths.vtp";
                int interval = 1;
                if (commands[i].parameters.count("interval") == 0)
                    interval = 1;
                else
                    interval = atoi(commands[i].parameters["interval"].c_str());
                trajs_vtk_pdt_to_vtp_3d(pathout+commands[i].parameters["filename"], atof(commands[i].parameters["color"].c_str()),interval);
            }

            if (commands[i].command == "initialize_marginal_v_dist")
            {
                show_in_window("initializing " +commands[i].parameters["dist"] + " distribution ..." );
                dist.name = commands[i].parameters["dist"];
                if (dist.name=="nonparameteric")
                {
                    show_in_window("reading distribution");
                    dist.readfromfile(pathout + commands[i].parameters["filename"]);
                }
                else
                {   for (int j = 0; j < split(commands[i].parameters["params"],',').size(); j++)
                        dist.params.push_back(atof(split(commands[i].parameters["params"],',')[j].c_str()));
                }

            }

            if (commands[i].command == "set_copula_params")
            {
                show_in_window("initializing copula..." );
                Copula.copula = commands[i].parameters["copula"];
                for (int j = 0; j < split(commands[i].parameters["params"],',').size(); j++)
                    Copula.parameters.push_back(atof(split(commands[i].parameters["params"],',')[j].c_str()));
                if (tolower(Copula.copula)=="gaussian")
                    Copula.SetCorrelation(Copula.parameters[0]);
                if (tolower(Copula.copula)=="frank")
                    Copula.Frank_copula_alpha = Copula.parameters[0];
                if (tolower(Copula.copula)=="experimental")
                {   TDMap copulamap;
                    copulamap.readfromfile(pathout + commands[i].parameters["experimental_copula_filename"]);
                    Copula.SetCopulaMap(copulamap);
                }


                Copula.SetDiffusionParams(atof(commands[i].parameters["diffusion"].c_str()),atof(commands[i].parameters["corr_ls"].c_str()),atof(commands[i].parameters["diffusion_corr_ls"].c_str()));

            }

            if (commands[i].command == "set_copula_diffusion_params")
            {
                show_in_window("initializing diffusion copula..." );
                Copula_diffusion.copula = commands[i].parameters["copula"];
                for (int j = 0; j < split(commands[i].parameters["params"],',').size(); j++)
                    Copula_diffusion.parameters.push_back(atof(split(commands[i].parameters["params"],',')[j].c_str()));
                if (tolower(Copula_diffusion.copula)=="gaussian")
                    Copula_diffusion.SetCorrelation(Copula_diffusion.parameters[0]);
                if (tolower(Copula_diffusion.copula)=="frank")
                    Copula_diffusion.Frank_copula_alpha = Copula_diffusion.parameters[0];
                if (tolower(Copula_diffusion.copula)=="experimental")
                {
                    TDMap copulamap;
                    copulamap.readfromfile(pathout + commands[i].parameters["experimental_copula_filename"]);
                    Copula_diffusion.SetCopulaMap(copulamap);
                }


            }


            if (commands[i].command == "create_correlated_v_pathways_ou")
            {
                show_in_window("Creating pathways ...");
                pset.create_ou_paths(atoi(commands[i].parameters["n"].c_str()), &dist, atof(commands[i].parameters["x_min"].c_str()), atof(commands[i].parameters["x_max"].c_str()), atof(commands[i].parameters["kappa"].c_str()), atof(commands[i].parameters["dx"].c_str()));

            }

            if (commands[i].command == "create_correlated_v_pathways_copula")
            {
                show_in_window("Creating pathways ...");
                pset.create_copula_paths(atoi(commands[i].parameters["n"].c_str()), &dist, atof(commands[i].parameters["x_min"].c_str()), atof(commands[i].parameters["x_max"].c_str()), atof(commands[i].parameters["epsilon"].c_str()), &Copula, atof(commands[i].parameters["dx"].c_str()));

            }


            if (commands[i].command == "write_pathways")
            {
                show_in_window("Writing pathways ...");
                pset.write(pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "write_pathways_vtp")
            {
                show_in_window("Writing pathways to vtp ...");
                vtkSmartPointer<vtkPolyDataMapper> mapper;
                mapper = pset.pathways_vtk_pdt_vtp(atof(commands[i].parameters["z_factor"].c_str()), atof(commands[i].parameters["offset"].c_str()));
                pset.write_vtk(mapper, pathout+commands[i].parameters["filename"]);
            }

            if (commands[i].command == "get_dist_at_time")
            {
                show_in_window("Getting distributions at time");
                dist_stores = pset.snapshotattime(atof(commands[i].parameters["t"].c_str()));
            }

            if (commands[i].command == "get_dist_at_location")
            {
                show_in_window("Getting distributions at location");
                dist_stores = pset.snapshotatlocation(atof(commands[i].parameters["x"].c_str()));
            }

            if (commands[i].command == "write_dist")
            {
                show_in_window("Writing distribution ...");
                dist_stores.write(pathout+commands[i].parameters["filename"].c_str());
            }

            if (commands[i].command == "write_dist_function")
            {
                show_in_window("writing distribution curve ...");
                if (commands[i].parameters.count("var") > 0)
                {
                    CBTC dist_curve = dist_stores.get_distribution(commands[i].parameters["var"], atoi(commands[i].parameters["log"].c_str()), atoi(commands[i].parameters["nbins"].c_str()));
                    dist_curve.writefile(pathout+commands[i].parameters["filename"]);
                    if (commands[i].parameters.count("filename_cum")>0)
                        dist_curve.getcummulative().writefile(pathout+commands[i].parameters["filename_cum"]);
                }
                else
                {
                    CBTCSet dist_curve = dist_stores.get_distribution(atoi(commands[i].parameters["log"].c_str()), atoi(commands[i].parameters["nbins"].c_str()));
                    dist_curve.writetofile(pathout+commands[i].parameters["filename"]);
                    if (commands[i].parameters.count("filename_cum")>0)
                        dist_curve.getcummulative().writetofile(pathout+commands[i].parameters["filename_cum"]);
                }


            }

            if (commands[i].command == "create_grid")
            {
                creategrid(atoi(commands[i].parameters["nx"].c_str()), atoi(commands[i].parameters["ny"].c_str()), atof(commands[i].parameters["dx"].c_str()), atof(commands[i].parameters["dy"].c_str()));
            }

            if (commands[i].command == "set_marginal_k_dist")
            {

                marginal_K_dist_type = commands[i].parameters["dist"];
                for (int j = 0; j < split(commands[i].parameters["params"],',').size(); j++)
                        marginal_K_dist_params.push_back(atof(split(commands[i].parameters["params"], ',')[j].c_str()));

            }

            if (commands[i].command == "make_paths_uniform")
            {
                show_in_window("Making trajectories uniform ... ");
                Traj.make_uniform_at_x(atof(commands[i].parameters["dx"].c_str()));
            }

            if (commands[i].command == "extract_pairs_diffusion")
            {

                if (!commands[i].parameters.count("nsequence"))
                    commands[i].parameters["nsequence"] = "2";

                CBTCSet pairs;
                double dt;
                bool fixed_interval = false;
                if (commands[i].parameters.count("fixed_interval")>0)
                    if (commands[i].parameters["fixed_interval"]=="true")
                        fixed_interval = true;


                dt = atof(commands[i].parameters["delta_t"].c_str());
                double D = atof(commands[i].parameters["diffusion"].c_str());
                pairs = get_correlation_based_on_random_samples_diffusion(atoi(commands[i].parameters["n"].c_str()),dt, D, dt/10, fixed_interval, commands[i].parameters["direction"]);
                double u_grad;

                pairs.writetofile(pathout+commands[i].parameters["filename"]);
                if (commands[i].parameters.count("dist_filename") > 0)
                {
                        pairs.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()),(pairs.BTC[0].maxC()-pairs.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).writefile(pathout+commands[i].parameters["dist_filename"]);
                }

                if (commands[i].parameters.count("v_grad_dist_file"))
                {
                    pairs.GetGradientDistribution(atof(commands[i].parameters["delta_t"].c_str()), atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["v_grad_dist_file"]);
                }

                if (commands[i].parameters.count("v_gnu_file"))
                {
                    TDMap GNU_out = pairs.get2DMap(atoi(commands[i].parameters["nbins"].c_str()));
                    GNU_out.writetofile_GNU(pathout + commands[i].parameters["v_gnu_file"],"", "v", "v'", "p(v,v')", true);
                }
                if (commands[i].parameters.count("ranks_filename") > 0)
                {

                    show_in_window("Writing ranks");
                    CBTCSet ranks(atoi(commands[i].parameters["nsequence"].c_str()));
                    bool _log = true;
                    if (commands[i].parameters.count("negative") > 0)
                        if (atoi(commands[i].parameters["negative"].c_str())==1)
                            _log = false;
                    show_in_window("negative = " + numbertostring(_log));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                    {
                        cout<<pairs.BTC[ii].n<<endl;
                        pairs.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),_log).writefile(pathout+commands[i].parameters["ranks_filename"]+"_"+numbertostring(ii));
                        cout<<pairs.BTC[ii].n<<endl;
                        ranks.BTC[ii] = pairs.BTC[ii].rank_bd(atoi(commands[i].parameters["nbins"].c_str()));
                        cout<<ranks.BTC[ii].n<<endl;
                        ranks.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["ranks_filename"]+"_u"+numbertostring(ii));
                    }
                    ranks.writetofile(pathout+commands[i].parameters["ranks_filename"]);
                    double alpha;
                    if (commands[i].parameters.count("copula_type")>0)
                        {

                            if (commands[i].parameters["copula_type"]=="frank")
                            {
                                show_in_window("Extracting Frank Copula Parameters");
                                CVector Likelihoods(2);
                                Likelihoods[0] = ranks.FrankCopulaLogLikelihood(2);
                                Likelihoods[1] = ranks.FrankCopulaLogLikelihood_deriv(2);
                                alpha = ranks.Estimate_Frank_Alpha();
                                extracted_OU_parameters.append("alpha",atof(commands[i].parameters["delta_t"].c_str()), alpha);
                                extracted_OU_parameters.append("correlation",atof(commands[i].parameters["delta_t"].c_str()), ranks.get_correlation());
                            }
                        }

                    if (commands[i].parameters.count("u_gnu_file"))
                    {
                        TDMap GNU_out = ranks.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["u_gnu_file"],"", "u", "u'", "p(u,u')");
                    }
                    if (commands[i].parameters.count("map_file"))
                    {
                        TDMap GNU_out = ranks.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile(pathout + commands[i].parameters["map_file"]);

                        if (commands[i].parameters.count("marginal_x_file_name")>0)
                        {
                            GNU_out.getcumulative("sym").writetofile(pathout+commands[i].parameters["marginal_x_file_name"]);
                        }

                        if (commands[i].parameters.count("theoretical_copula_filename")>0)
                        {
                            if (commands[i].parameters.count("copula_type")>0)
                            {
                                if (commands[i].parameters["copula_type"]=="frank")
                                {
                                    show_in_window("writing theoretical Frank copula density");
                                    CCopula FrankCopula;
                                    FrankCopula.copula = "frank";
                                    FrankCopula.Frank_copula_alpha = alpha;
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&FrankCopula);
                                }
                                if (commands[i].parameters["copula_type"]=="experimental")
                                {
                                    show_in_window("writing theoretical Experimental copula density");
                                    CCopula ExperimentalCopula;
                                    ExperimentalCopula.copula = "experimental";
                                    ExperimentalCopula.SetCopulaMap(GNU_out);
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&ExperimentalCopula);
                                    if (commands[i].parameters.count("as_array_filename")>0)
                                    {
                                        GNU_out.writetheoreticalcopulatofile_points(pathout + commands[i].parameters["as_array_filename"],&ExperimentalCopula);
                                    }
                                }
                            }
                        }

                    }
                    if (commands[i].parameters.count("joint_cdf_file"))
                    {
                        TDMap GNU_out = ranks.getJointCDF(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile(pathout + commands[i].parameters["joint_cdf_file"]);
                    }
                    u_grad = diff(ranks.BTC[0],ranks.BTC[1])/pow(atof(commands[i].parameters["delta_t"].c_str()),2)/double(ranks.BTC[0].n);
                }
                if (commands[i].parameters.count("normal_filename") > 0)
                {
                    show_in_window("Writing normals");
                    CBTCSet normals(atoi(commands[i].parameters["nsequence"].c_str()));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                        normals.BTC[ii] = pairs.BTC[ii].map_to_standard_normal(atoi(commands[i].parameters["nbins"].c_str()));

                    if (commands[i].parameters.count("w_gnu_file"))
                    {
                        TDMap GNU_out = normals.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),-4,4);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["w_gnu_file"],"", "w", "w'", "p(w,w')");
                    }

                    if (commands[i].parameters.count("w_grad_dist_file"))
                    {
                        normals.GetGradientDistribution(atof(commands[i].parameters["delta_t"].c_str()),atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["w_grad_dist_file"]);
                    }

                    normals.writetofile(pathout + commands[i].parameters["normal_filename"]);
                    if (commands[i].parameters.count("OU_parameters_filename") > 0)
                    {
                        bool gaussian_copula = false;
                            if (commands[i].parameters.count("copula_type")>0)
                                if (commands[i].parameters["copula_type"]!="gaussian")
                                    gaussian_copula = false;
                        if (gaussian_copula)
                        {    show_in_window("Calculating OU params");


                            double corr = normals.get_correlation();
                            double standard_deviation = (normals.BTC[0]-normals.BTC[1]).moment(2);

                            extracted_OU_parameters.append("Ro", atof(commands[i].parameters["delta_t"].c_str()), corr);
                            extracted_OU_parameters.append("u_grad", atof(commands[i].parameters["delta_t"].c_str()),u_grad);
                            extracted_OU_parameters.append("moment2", atof(commands[i].parameters["delta_t"].c_str()),standard_deviation);
                            CVector X(3);
                            X[0] = corr;
                            X[1] = u_grad;
                            X[2] = standard_deviation;
                            show_in_window("Writing OU params");
                            X.writetofile(pathout + commands[i].parameters["OU_parameters_filename"]);
                            if (commands[i].parameters.count("theoretical_copula_filename")>0)
                            {   if (gaussian_copula)
                                {
                                    TDMap GNU_out(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                                    show_in_window("writing theoretical Gaussian copula density");
                                    CCopula GaussianCopula;
                                    GaussianCopula.copula = "gaussian";
                                    GaussianCopula.SetCorrelation(-log(corr));
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&GaussianCopula);
                                }
                            }

                        }


                    }
                }
            }

            if (commands[i].command == "readtrajsfromfile")
            {
                show_in_window("Reading trajectories from file '" + commands[i].parameters["filename"] + "'");
                if (!Traj.getfromMODflowfile(commands[i].parameters["filename"]))
                {
                    show_in_window("Reading trajectories failed");
                }
                else
                    show_in_window("Reading trajectories completed");
            }

            if (commands[i].command == "readtrajsfromfracfile")
            {
                show_in_window("Reading trajectories from file '" + commands[i].parameters["filename"] + "'");
                int column_number = 0;
                string reactionfile = "";
                if (commands[i].parameters.count("reaction_state_file")>0)
                {
                    reactionfile = commands[i].parameters["reaction_state_file"];
                }

                if (commands[i].parameters.count("column")>0)
                {
                    column_number = atoi(commands[i].parameters["column"].c_str());
                }

                if (!Traj.getfromShermanfile(commands[i].parameters["filename"],reactionfile, column_number))
                {
                    show_in_window("Reading trajectories failed");
                }
                else
                    show_in_window("Reading trajectories completed");
            }

            if (commands[i].command == "readtrajsfromfracfile_v")
            {
                show_in_window("Reading trajectories from file '" + commands[i].parameters["filename"] + "'");
                if (!Traj.getfromShermanfile_v(commands[i].parameters["filename"]))
                {
                    show_in_window("Reading trajectories failed");
                }
                else
                    show_in_window("Reading trajectories completed");
            }

            if (commands[i].command == "extract_pairs")
            {
                show_in_window("Extracting pairs ... ");
                CBTCSet pairs;
                bool random_sampling = false;
                bool time_based = false;

                if (commands[i].parameters.count("random_sampling")==1)
                    if (atoi(commands[i].parameters["random_sampling"].c_str())==1)
                        random_sampling = true;

                if (commands[i].parameters.count("time_based")==1)
                    if (atoi(commands[i].parameters["time_based"].c_str())==1)
                        time_based = true;


                if (!commands[i].parameters.count("nsequence"))
                            commands[i].parameters["nsequence"] = "2";

                bool magnitude = false;
                if (commands[i].parameters.count("magnitude"))
                {
                    magnitude = atoi(commands[i].parameters["magnitude"].c_str());
                }

                cout<<"Magnitude = "<<magnitude<<endl;
                double dx;
                if (random_sampling)
                    {
                        if (!time_based)
                        {

                            dx = atof(commands[i].parameters["delta_x"].c_str());
                            pairs = get_correlation_based_on_random_samples(atoi(commands[i].parameters["n"].c_str()),dx,atof(commands[i].parameters["increment"].c_str()),magnitude);
                        }
                        else
                        {
                            double dt;
                            dx = dt = atof(commands[i].parameters["delta_t"].c_str());
                            pairs = get_correlation_based_on_random_samples_dt(atoi(commands[i].parameters["n"].c_str()),dt,atof(commands[i].parameters["increment"].c_str()));
                        }

                    }
                    else
                    {

                        if (Traj.paths.size())
                            pairs = Traj.get_pair_v(atoi(commands[i].parameters["increment"].c_str()), atoi(commands[i].parameters["n"].c_str()),atoi(commands[i].parameters["nsequence"].c_str()));
                        else if (pset.paths.size())
                            pairs = pset.get_pair_v(atoi(commands[i].parameters["increment"].c_str()), atoi(commands[i].parameters["n"].c_str()),atoi(commands[i].parameters["nsequence"].c_str()));
                        else
                        {
                            show_in_window("No trajectories has been initialized!");
                            return;
                        }
                    }



                pairs.writetofile(pathout+commands[i].parameters["filename"]);
                if (commands[i].parameters.count("dist_filename") > 0)
                {
                    pairs.BTC[0].distribution(atoi(commands[i].parameters["nbins"].c_str()),(pairs.BTC[0].maxC()-pairs.BTC[0].minC())*atof(commands[i].parameters["smoothing_factor"].c_str())).writefile(pathout+commands[i].parameters["dist_filename"]);
                }

                if (commands[i].parameters.count("v_grad_dist_file"))
                {
                    pairs.GetGradientDistribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["v_grad_dist_file"]);
                }

                if (commands[i].parameters.count("v_gnu_file"))
                {
                    TDMap GNU_out = pairs.get2DMap(atoi(commands[i].parameters["nbins"].c_str()));
                    GNU_out.writetofile_GNU(pathout + commands[i].parameters["v_gnu_file"],"", "v", "v'", "p(v,v')", true);
                }

                bool _log = true;


                if (commands[i].parameters.count("ranks_filename") > 0)
                {
                    show_in_window("Writing ranks");
                    CBTCSet ranks(atoi(commands[i].parameters["nsequence"].c_str()));
                      if (commands[i].parameters.count("negative") > 0)
                        if (atoi(commands[i].parameters["negative"].c_str())==1)
                            _log = false;
                    show_in_window("negative = " + numbertostring(_log));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                    {
                        cout<<pairs.BTC[ii].n<<endl;
                        pairs.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),_log).writefile(pathout+commands[i].parameters["ranks_filename"]+"_"+numbertostring(ii));
                        cout<<pairs.BTC[ii].n<<endl;
                        ranks.BTC[ii] = pairs.BTC[ii].rank_bd(atoi(commands[i].parameters["nbins"].c_str()),_log);
                        cout<<ranks.BTC[ii].n<<endl;
                        ranks.BTC[ii].getcummulative_direct(atoi(commands[i].parameters["nbins"].c_str()),_log).writefile(pathout+commands[i].parameters["ranks_filename"]+"_u"+numbertostring(ii));
                    }
                    ranks.writetofile(pathout+commands[i].parameters["ranks_filename"]);
                    double alpha;
                    if (commands[i].parameters.count("copula_type")>0)
                        {

                            if (commands[i].parameters["copula_type"]=="frank")
                            {
                                show_in_window("Extracting Frank Copula Parameters");
                                CVector Likelihoods(2);
                                alpha = ranks.Estimate_Frank_Alpha();
                                extracted_OU_parameters.append("alpha",dx, alpha);
                                extracted_OU_parameters.append("correlation",atof(commands[i].parameters["delta_t"].c_str()), ranks.get_correlation());
                            }
                        }
                    if (commands[i].parameters.count("u_gnu_file"))
                    {
                        TDMap GNU_out = ranks.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["u_gnu_file"],"", "u", "u'", "p(u,u')");
                    }
                    if (commands[i].parameters.count("map_file"))
                    {
                        TDMap GNU_out = ranks.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile(pathout + commands[i].parameters["map_file"]);

                        if (commands[i].parameters.count("theoretical_copula_filename")>0)
                        {
                            if (commands[i].parameters.count("copula_type")>0)
                            {
                                if (commands[i].parameters["copula_type"]=="frank")
                                {
                                    show_in_window("writing theoretical Frank copula density");
                                    CCopula FrankCopula;
                                    FrankCopula.copula = "frank";
                                    FrankCopula.Frank_copula_alpha = alpha;
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&FrankCopula);
                                }
                                if (commands[i].parameters["copula_type"]=="experimental")
                                {
                                    show_in_window("writing theoretical experimental copula density");
                                    CCopula ExperimentalCopula;
                                    ExperimentalCopula.copula = "experimental";
                                    ExperimentalCopula.SetCopulaMap(GNU_out);
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&ExperimentalCopula);
                                    if (commands[i].parameters.count("as_array_filename")>0)
                                    {
                                        GNU_out.writetheoreticalcopulatofile_points(pathout + commands[i].parameters["as_array_filename"],&ExperimentalCopula);
                                    }
                                }
                                if (commands[i].parameters["copula_type"]=="gaussian")
                                {
                                    show_in_window("writing theoretical experimental copula density");
                                    CCopula ExperimentalCopula;
                                    ExperimentalCopula.copula = "gaussian";
                                    ExperimentalCopula.SetCopulaMap(GNU_out);
                                    GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&ExperimentalCopula);
                                    if (commands[i].parameters.count("as_array_filename")>0)
                                    {
                                        GNU_out.writetheoreticalcopulatofile_points(pathout + commands[i].parameters["as_array_filename"],&ExperimentalCopula);
                                    }
                                }
                            }
                        }
                    }


                    if (commands[i].parameters.count("joint_cdf_file"))
                    {
                        TDMap GNU_out = ranks.getJointCDF(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                        GNU_out.writetofile(pathout + commands[i].parameters["joint_cdf_file"]);
                    }
                }
                if (commands[i].parameters.count("normal_filename") > 0)
                {
                    show_in_window("Writing normals");
                    CBTCSet normals(atoi(commands[i].parameters["nsequence"].c_str()));
                    for (int ii=0; ii<atoi(commands[i].parameters["nsequence"].c_str()); ii++)
                        normals.BTC[ii] = pairs.BTC[ii].map_to_standard_normal(atoi(commands[i].parameters["nbins"].c_str()));

                    if (commands[i].parameters.count("w_gnu_file"))
                    {
                        TDMap GNU_out = normals.get2DMap(atoi(commands[i].parameters["nbins"].c_str()),-4,4);
                        GNU_out.writetofile_GNU(pathout + commands[i].parameters["w_gnu_file"],"", "w", "w'", "p(w,w')");
                    }


                    if (commands[i].parameters.count("w_grad_dist_file"))
                    {
                        normals.GetGradientDistribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["w_grad_dist_file"]);
                    }


                    normals.writetofile(pathout + commands[i].parameters["normal_filename"]);
                    bool gaussian_copula = true;
                        if (commands[i].parameters.count("copula_type")>0)
                        {    if (commands[i].parameters["copula_type"]!="gaussian")
                                gaussian_copula = false;
                        }
                    if (commands[i].parameters.count("OU_parameters_filename") > 0)
                    {
                        if (gaussian_copula)
                        {
                            show_in_window("Calculating OU params");

                            if (!time_based)
                            {

                                show_in_window("Writing OU parameters...");
                                CVector X = normals.get_kappa_gamma(atof(commands[i].parameters["delta_x"].c_str()));

                                extracted_OU_parameters.append("Gamma", atof(commands[i].parameters["delta_x"].c_str()), X[0]);
                                extracted_OU_parameters.append("Kappa", atof(commands[i].parameters["delta_x"].c_str()), X[1]);
                                double corr = normals.get_correlation();
                                extracted_OU_parameters.append("Ro", atof(commands[i].parameters["delta_x"].c_str()), corr);
                                if (atoi(commands[i].parameters["nsequence"].c_str())>2)
                                    extracted_OU_parameters.append("p3", atoi(commands[i].parameters["increment"].c_str()), X[2]);

                                show_in_window("Writing OU params");
                                X.writetofile(pathout + commands[i].parameters["OU_parameters_filename"]);

                            if (commands[i].parameters.count("theoretical_copula_filename")>0)
                                {   if (gaussian_copula)
                                    {
                                        TDMap GNU_out(atoi(commands[i].parameters["nbins"].c_str()),0,1);
                                        show_in_window("writing theoretical Gaussian copula density");
                                        CCopula GaussianCopula;
                                        GaussianCopula.copula = "gaussian";
                                        GaussianCopula.SetCorrelation(-log(corr));
                                        GNU_out.writetheoreticalcopulatofile(pathout + commands[i].parameters["theoretical_copula_filename"],&GaussianCopula);
                                    }
                                }
                            }
                            else
                            {
                                CVector X = normals.get_kappa_gamma(atof(commands[i].parameters["delta_t"].c_str()));
                                extracted_OU_parameters.append("Gamma", atof(commands[i].parameters["delta_t"].c_str()), X[0]);
                                extracted_OU_parameters.append("Kappa", atof(commands[i].parameters["delta_t"].c_str()), X[1]);

                                show_in_window("Writing OU params");
                                X.writetofile(pathout + commands[i].parameters["OU_parameters_filename"]);
                            }
                        }
                    }
                }
            }

            if (commands[i].command == "write_all_ou_parameters")
            {
                show_in_window("Writing OU parameters ... ");
                extracted_OU_parameters.writetofile(pathout+commands[i].parameters["filename"]);
                if (commands[i].parameters.count("density_filename") > 0)
                {
                        CBTCSet OU_Density;
                        for (int j = 0; j < extracted_OU_parameters.nvars; j++)
                        {   CBTC BTC = extracted_OU_parameters.BTC[j].distribution(atoi(commands[i].parameters["nbins"].c_str()),0);
                            OU_Density.append(BTC,extracted_OU_parameters.BTC[j].name);

                        }

                        OU_Density.writetofile(pathout+commands[i].parameters["density_filename"]);
                }
            }

            if (commands[i].command == "clear_all")
            {
                show_in_window("Clearing variables ... ");
                for (int i=0; i<C.size(); i++)
                    C[i].matr.clear();
                H.matr.clear();
                KD.matr.clear();
                Kt.matr.clear();
                Kv.matr.clear();
                clear();
                trajectories.clear();
                Traj.paths.clear();
                sect_dist.clear();
                pts.clear();
            }

            if (commands[i].command == "write")
            {
                    show_in_window(commands[i].parameters["content"]);
            }
            if (commands[i].command == "read_3d_grid_data")
            {
                threedgrid.ReadFromModFlowFile(commands[i].parameters["filename"]);
            }
            if (commands[i].command == "write_3d_grid_to_vts")
            {
                threedgrid.create_vts(pathout + commands[i].parameters["filename"], atoi(commands[i].parameters["x_limit"].c_str()), atoi(commands[i].parameters["y_limit"].c_str()), atoi(commands[i].parameters["z_limit"].c_str()),atoi(commands[i].parameters["x_interval"].c_str()), atoi(commands[i].parameters["y_interval"].c_str()), atoi(commands[i].parameters["z_interval"].c_str()));
            }
            if (commands[i].command == "get_residuals")
            {
                CTimeSeries out = GetResiduals(commands[i].parameters["property"]);
                if (commands[i].parameters.count("raw_filename")>0)
                {
                    out.writefile(pathout+commands[i].parameters["raw_filename"]);
                }
                if (commands[i].parameters.count("dist_filename")>0)
                {
                    out.distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["dist_filename"]);
                }
                if (commands[i].parameters.count("std_filename")>0)
                {
                    CVector STD(2);
                    STD[0] = GP.dx;
                    STD[1] = out.std();
                    STD.writetofile(pathout+commands[i].parameters["std_filename"]);
                }
            }
            if (commands[i].command == "get_v_residuals")
            {
                CTimeSeries out = GetResiduals("omega");
                if (commands[i].parameters.count("raw_filename")>0)
                {
                    out.writefile(pathout+commands[i].parameters["raw_filename"]);
                }
                if (commands[i].parameters.count("dist_filename")>0)
                {
                    out.distribution(atoi(commands[i].parameters["nbins"].c_str())).writefile(pathout+commands[i].parameters["dist_filename"]);
                }
                if (commands[i].parameters.count("std_filename")>0)
                {
                    CVector STD(2);
                    STD[0] = GP.dx;
                    STD[1] = out.std();
                    STD.writetofile(pathout+commands[i].parameters["std_filename"]);
                }
            }

	}

	show_in_window("Finished!");
#ifdef QT_version
	main_window->setCursor(Qt::ArrowCursor);
#endif // QT_version
}
